import math
import json
import pandas as pd
import numpy as np
import networkx as nx
import time
from ortools.sat.python import cp_model

# ==========================================
# 1. FUNCIONES AUXILIARES
# ==========================================

def haversine_distance(lon1, lat1, lon2, lat2):
    """Calcula distancia real en km entre dos puntos geográficos."""
    R = 6371.0 
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

# ==========================================
# 2. PRE-PROCESAMIENTO
# ==========================================

def preprocess_data(csv_path, geojson_path, pruning_budget_k=None):
    print(f"⚙️ [1/2] Pre-procesando datos (Grafo + Dijkstra)...")
    start_time = time.time()
    
    # Cargar datos
    df = pd.read_csv(csv_path)
    
    SCALE_COST = 1000
    SCALE_SCORE = 10
    SCALE_AREA = 100
    
    ids = df['grid_id'].values
    n_nodes = len(ids)
    id_map = {original_id: i for i, original_id in enumerate(ids)}
    
    S_LIST = ['atelerix', 'martes', 'eliomys', 'oryctolagus']
    n_species = len(S_LIST)
    W_vals = [1.0, 1.0, 2.0, 1.5]
    
    areas = (df['cell_area_km2'].values * SCALE_AREA).astype(int).tolist()
    suitability = [df[f'suitability_{s}'].values.tolist() for s in S_LIST]
    cost_adapt = [(df[f'cost_adaptation_{s}'].values * SCALE_COST).astype(int).tolist() for s in S_LIST]
    
    col_map = {'atelerix': 'has_atelerix_algirus', 'martes': 'has_martes_martes',
               'eliomys': 'has_eliomys_quercinus', 'oryctolagus': 'has_oryctolagus_cuniculus'}
    has_spec = [df[col_map[s]].values.astype(bool).tolist() for s in S_LIST]

    # Coordenadas Reales
    with open(geojson_path, 'r') as f:
        geo_data = json.load(f)
    centroids = {}
    for feature in geo_data['features']:
        gid = feature['properties']['grid_id']
        geom = feature['geometry']
        if geom['type'] == 'Polygon':
            coords = geom['coordinates'][0]
        elif geom['type'] == 'MultiPolygon':
            coords = geom['coordinates'][0][0]
        lon = sum(p[0] for p in coords)/len(coords)
        lat = sum(p[1] for p in coords)/len(coords)
        centroids[gid] = (lon, lat)
    
    lons = np.array([centroids.get(i, (0,0))[0] for i in ids])
    lats = np.array([centroids.get(i, (0,0))[1] for i in ids])

    # Construir Grafo
    G = nx.Graph()
    cost_corr_base = df['cost_corridor'].values
    neighbors_series = df['neighbors'].astype(str).str.split(r'[;,]')
    
    edges_unique_cost = {} 
    all_edge_costs = []

    for u_idx, n_list in enumerate(neighbors_series):
        if not isinstance(n_list, list): continue
        for v_str in n_list:
            if v_str in id_map:
                v_idx = id_map[v_str]
                if u_idx < v_idx: 
                    dist_km = haversine_distance(lons[u_idx], lats[u_idx], lons[v_idx], lats[v_idx])
                    avg_friction = (cost_corr_base[u_idx] + cost_corr_base[v_idx]) / 2
                    cost = int(avg_friction * dist_km * SCALE_COST)
                    
                    G.add_edge(u_idx, v_idx, weight=cost)
                    edges_unique_cost[(u_idx, v_idx)] = cost
                    all_edge_costs.append(cost)

    # Límite de Poda
    cutoff_value = None
    if pruning_budget_k:
        median_cost = np.median(all_edge_costs) if all_edge_costs else 1
        budget_int = int(pruning_budget_k * SCALE_COST)
        limit_pct = budget_int * 0.15 
        safety_floor = median_cost * 6
        cutoff_value = max(limit_pct, safety_floor)

    # Dijkstra Multi-Source
    paths_requirements = [[None for _ in range(n_species)] for _ in range(n_nodes)]
    
    for s_idx in range(n_species):
        sources = [i for i, val in enumerate(has_spec[s_idx]) if val]
        if not sources: continue
        
        dists, paths = nx.multi_source_dijkstra(G, sources, weight='weight', cutoff=cutoff_value)

        for target, path_nodes in paths.items():
            if target in sources:
                paths_requirements[target][s_idx] = []
                continue
            edge_list = []
            for k in range(len(path_nodes) - 1):
                u, v = sorted((path_nodes[k], path_nodes[k+1]))
                edge_list.append((u, v))
            paths_requirements[target][s_idx] = edge_list

    elapsed = time.time() - start_time
    print(f"✅ Pre-procesamiento completado en {elapsed:.2f} s.")
    
    return {
        'ids': ids, 'n_nodes': n_nodes, 'S_LIST': S_LIST, 'W_vals': W_vals,
        'areas': areas, 'suitability': suitability, 'cost_adapt': cost_adapt, 'has_spec': has_spec,
        'edges_unique_cost': edges_unique_cost, 'paths_requirements': paths_requirements,
        'SCALE_COST': SCALE_COST, 'SCALE_SCORE': SCALE_SCORE, 'SCALE_AREA': SCALE_AREA
    }

# ==========================================
# 3. SOLVER (CP-SAT)
# ==========================================

def solve_conservation(data, budget=500, time_limit_sec=600, n_threads=12):
    print(f"\n⚡ [2/2] INICIANDO SOLVER (Presupuesto: {budget} k€)")
    solver_start_time = time.time()
    
    # Desempaquetar
    n_nodes = data['n_nodes']; S_LIST = data['S_LIST']; n_species = len(S_LIST)
    W_vals = data['W_vals']; areas = data['areas']
    suitability = data['suitability']; cost_adapt = data['cost_adapt']
    has_spec = data['has_spec']; edges_unique_cost = data['edges_unique_cost']
    paths_requirements = data['paths_requirements']
    
    SCALE_COST = data['SCALE_COST']
    SCALE_SCORE = data['SCALE_SCORE']
    SCALE_AREA = data['SCALE_AREA']
    budget_int = int(budget * SCALE_COST)

    model = cp_model.CpModel()

    # Variables
    x = [[model.NewBoolVar(f'x_{i}_{s}') for s in range(n_species)] for i in range(n_nodes)]
    y = [[model.NewBoolVar(f'y_{i}_{s}') for s in range(n_species)] for i in range(n_nodes)]
    stress = [model.NewBoolVar(f'strs_{i}') for i in range(n_nodes)]

    # Z: Variables de arista
    used_edges = set()
    for i in range(n_nodes):
        for s in range(n_species):
            reqs = paths_requirements[i][s]
            if reqs: used_edges.update(reqs)
            
    z = {}
    for e in used_edges:
        z[e] = model.NewBoolVar(f'z_{e}')

    # Restricciones
    obj_terms = []
    cost_terms = []
    active_area_terms = [[] for _ in range(n_species)]
    PENALTY_STRESS = int(350.0 * SCALE_SCORE)

    for i in range(n_nodes):
        model.Add(x[i][1] + x[i][2] <= 1)
        model.Add(stress[i] >= x[i][1] + x[i][3] - 1)
        obj_terms.append(stress[i] * (-PENALTY_STRESS))

        for s in range(n_species):
            req_edges = paths_requirements[i][s]
            
            # Poda
            if req_edges is None:
                model.Add(x[i][s] == 0)
                continue

            model.Add(y[i][s] <= x[i][s])
            if not has_spec[s][i]: model.Add(x[i][s] <= y[i][s])
            cost_terms.append(y[i][s] * cost_adapt[s][i])

            if req_edges:
                for e in req_edges:
                    model.AddImplication(x[i][s], z[e])

            coef_x = int((suitability[s][i] * W_vals[s] * areas[i] * SCALE_SCORE) / SCALE_AREA)
            obj_terms.append(x[i][s] * coef_x)
            
            added = 3.0 - suitability[s][i]
            if added > 0:
                coef_y = int((added * W_vals[s] * areas[i] * SCALE_SCORE) / SCALE_AREA)
                obj_terms.append(y[i][s] * (coef_y - 1))
            else:
                obj_terms.append(y[i][s] * -1)
                
            active_area_terms[s].append(x[i][s] * areas[i])

    for e, z_var in z.items():
        cost = edges_unique_cost.get(e, 0)
        cost_terms.append(z_var * cost)
        obj_terms.append(z_var * -1) 

    if budget < 100_000: model.Add(sum(cost_terms) <= budget_int)

    total_active = sum(sum(active_area_terms[s]) for s in range(n_species))
    min_pct = [5, 20, 5, 15]
    for s in range(n_species):
        s_sum = sum(active_area_terms[s])
        model.Add(s_sum * 100 >= min_pct[s] * total_active)

    model.Maximize(sum(obj_terms))

    # Ejecución
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = float(time_limit_sec)
    solver.parameters.num_search_workers = n_threads
    # solver.parameters.relative_gap_tolerance = 0.01 # Desactivado para buscar óptimo
    solver.parameters.log_search_progress = False
    
    print("🧠 Buscando óptimo...")
    status = solver.Solve(model)
    
    total_time = time.time() - solver_start_time

    if status in [cp_model.OPTIMAL, cp_model.FEASIBLE]:
        obj_val = solver.ObjectiveValue()
        best_bound = solver.BestObjectiveBound()
        gap = 0.0
        if best_bound != 0:
            gap = abs(best_bound - obj_val) / abs(best_bound) * 100
            
        # Preparar datos de salida
        res = []
        for i in range(n_nodes):
            row = {'grid_id': data['ids'][i]}
            for s_idx, s_name in enumerate(S_LIST):
                if paths_requirements[i][s_idx] is not None:
                    row[f'active_{s_name}'] = solver.Value(x[i][s_idx])
                    row[f'invest_{s_name}'] = solver.Value(y[i][s_idx])
                else:
                    row[f'active_{s_name}'] = 0
                    row[f'invest_{s_name}'] = 0
            res.append(row)
            
        df_sol = pd.DataFrame(res)
        
        # Diccionario de Métricas para el print final
        metrics = {
            'score': obj_val / SCALE_SCORE,
            'gap': gap,
            'time': total_time,
            'status': solver.StatusName(status)
        }
        return df_sol, metrics
    else:
        print(f"❌ Fallo: {solver.StatusName(status)}")
        return None, None

# ==========================================
# 4. EJECUCIÓN Y REGISTRO
# ==========================================
if __name__ == "__main__":
    CSV_FILE = 'final_dataset.csv'
    GEOJSON_FILE = 'final_dataset.geojson'
    PRESUPUESTO = 500
    N_RUNS = 10
    
    # 1. Pre-procesar (Solo una vez, es determinista)
    datos = preprocess_data(CSV_FILE, GEOJSON_FILE, pruning_budget_k=PRESUPUESTO)
    
    print(f"🚀 Iniciando análisis de robustez ({N_RUNS} ejecuciones)...")
    threads = [1] #Los que se quieran probar, e.g., [1, 4, 8, 12]

    for n_threads in threads:
        print(f"\n=== Usando {n_threads} hilos de búsqueda ===")

        # CAMBIO 2: Resetear resultados por configuración de hilos
        resultados_globales = pd.DataFrame()

        for i in range(N_RUNS):
            print(f"\n🔄 Ejecución {i+1}/{N_RUNS} ...")
            
            # Llamamos al solver con una semilla distinta cada vez
            df_sol, metrics = solve_conservation(
                datos, 
                budget=PRESUPUESTO, 
                time_limit_sec=1200, 
                n_threads=n_threads
            )
            
            if metrics is not None:
                metrics['run_id'] = i

                # 🔥 CAMBIO 1: Acumular correctamente una fila por ejecución
                fila = {m: metrics[m] for m in metrics.keys()}
                resultados_globales = pd.concat(
                    [resultados_globales, pd.DataFrame([fila])],
                    ignore_index=True
                )
                
                # Guardado de resultados por hilos
                filename = f'cpsat_{n_threads}_threads_{PRESUPUESTO}.csv'
                resultados_globales.to_csv(filename, index=False)

        # 3. Análisis de Resultados
        df_res = pd.DataFrame(resultados_globales)
        
        print("\n" + "="*40)
        print("📊 RESUMEN ESTADÍSTICO (30 RUNS)")
        print("="*40)
        print(df_res[['score', 'gap', 'time']].describe().round(2))