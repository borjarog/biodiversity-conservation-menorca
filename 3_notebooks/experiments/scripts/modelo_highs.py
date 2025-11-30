import math
import json
import pandas as pd
import numpy as np
import networkx as nx
import time
from ortools.linear_solver import pywraplp

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
# 2. PRE-PROCESAMIENTO (NetworkX + Poda)
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
    
    # Datos Numéricos
    areas = (df['cell_area_km2'].values * SCALE_AREA).astype(int).tolist()
    suitability = [df[f'suitability_{s}'].values.tolist() for s in S_LIST]
    cost_adapt = [(df[f'cost_adaptation_{s}'].values * SCALE_COST).astype(int).tolist() for s in S_LIST]
    
    col_map = {'atelerix': 'has_atelerix_algirus', 'martes': 'has_martes_martes',
               'eliomys': 'has_eliomys_quercinus', 'oryctolagus': 'has_oryctolagus_cuniculus'}
    has_spec = [df[col_map[s]].values.astype(bool).tolist() for s in S_LIST]

    # Coordenadas Reales (GeoJSON)
    with open(geojson_path, 'r') as f:
        geo_data = json.load(f)
    centroids = {}
    for feature in geo_data['features']:
        gid = feature['properties']['grid_id']
        geom = feature['geometry']
        if geom['type'] == 'Polygon': coords = geom['coordinates'][0]
        elif geom['type'] == 'MultiPolygon': coords = geom['coordinates'][0][0]
        lon = sum(p[0] for p in coords)/len(coords); lat = sum(p[1] for p in coords)/len(coords)
        centroids[gid] = (lon, lat)
    
    lons = np.array([centroids.get(i, (0,0))[0] for i in ids])
    lats = np.array([centroids.get(i, (0,0))[1] for i in ids])

    # Construir Grafo NetworkX
    G = nx.Graph()
    cost_corr_base = df['cost_corridor'].values
    neighbors_series = df['neighbors'].astype(str).str.split(r'[;,]')
    
    # --- IMPORTANTE: Esta variable 'adj' es la que faltaba en el return ---
    adj = [[] for _ in range(n_nodes)]
    
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
                    
                    # Guardamos la adyacencia para reconstruir corredores luego
                    adj[u_idx].append((v_idx, cost))
                    adj[v_idx].append((u_idx, cost))

    # Límite de Poda
    cutoff_value = None
    if pruning_budget_k:
        median_cost = np.median(all_edge_costs) if all_edge_costs else 1
        budget_int = int(pruning_budget_k * SCALE_COST)
        limit_pct = budget_int * 0.15 
        safety_floor = median_cost * 6
        cutoff_value = max(limit_pct, safety_floor)
        print(f"   ✂️ Poda Activada: Ignorando rutas > {cutoff_value/SCALE_COST:.2f} k€")

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
        # AÑADIDO 'adj' AQUÍ:
        'adj': adj,
        'edges_unique_cost': edges_unique_cost, 'paths_requirements': paths_requirements,
        'SCALE_COST': SCALE_COST, 'SCALE_SCORE': SCALE_SCORE, 'SCALE_AREA': SCALE_AREA
    }

# ==========================================
# 3. SOLVER (MIP - HIGHS)
# ==========================================
def solve_conservation_highs(data, budget=1000, time_limit_sec=600):
    print(f"\n⚡ [2/2] INICIANDO SOLVER HIGHS (Presupuesto: {budget} k€, 8 Cores)")
    solver_start_time = time.time()
    
    # 1. Crear Solver HiGHS
    try:
         solver = pywraplp.Solver.CreateSolver('HIGHS_MIP')
         if not solver: solver = pywraplp.Solver.CreateSolver('HIGHS')
    except:
         solver = pywraplp.Solver.CreateSolver('HIGHS')

    if not solver:
        print(f"❌ Error: El solver HiGHS no está disponible.")
        return None

    # 2. Configuración
    solver.SetTimeLimit(int(time_limit_sec * 1000)) # milisegundos
    solver.SetNumThreads(12) # REQUISITO: 8 NUCLEOS
    
    # Desempaquetar
    n_nodes = data['n_nodes']; S_LIST = data['S_LIST']; n_species = len(S_LIST)
    W_vals = data['W_vals']; areas = data['areas']
    suitability = data['suitability']; cost_adapt = data['cost_adapt']
    has_spec = data['has_spec']; edges_unique_cost = data['edges_unique_cost']
    paths_requirements = data['paths_requirements']
    
    SCALE_COST = data['SCALE_COST']; SCALE_SCORE = data['SCALE_SCORE']; SCALE_AREA = data['SCALE_AREA']
    budget_int = int(budget * SCALE_COST)

    # --- VARIABLES (IntVar 0-1 para MIP) ---
    x = {}; y = {}; stress = {}
    for i in range(n_nodes):
        stress[i] = solver.IntVar(0, 1, f"stress_{i}")
        for s in range(n_species):
            x[i, s] = solver.IntVar(0, 1, f"x_{i}_{s}")
            y[i, s] = solver.IntVar(0, 1, f"y_{i}_{s}")

    used_edges = set()
    for i in range(n_nodes):
        for s in range(n_species):
            reqs = paths_requirements[i][s]
            if reqs: used_edges.update(reqs)
            
    z = {}
    for e in used_edges:
        z[e] = solver.IntVar(0, 1, f"z_{e}")

    # --- RESTRICCIONES ---
    obj_expr = 0
    cost_expr = 0
    active_area_exprs = [0 for _ in range(n_species)]
    PENALTY_STRESS = int(350.0 * SCALE_SCORE)

    for i in range(n_nodes):
        # Stress
        solver.Add(x[i, 1] + x[i, 2] <= 1)
        solver.Add(stress[i] >= x[i, 1] + x[i, 3] - 1)
        obj_expr += stress[i] * (-PENALTY_STRESS)

        for s in range(n_species):
            req_edges = paths_requirements[i][s]
            
            # Poda
            if req_edges is None:
                solver.Add(x[i, s] == 0)
                continue

            # Inversión y <= x
            solver.Add(y[i, s] <= x[i, s])
            if not has_spec[s][i]: solver.Add(x[i, s] <= y[i, s])
            
            cost_expr += y[i, s] * cost_adapt[s][i]

            # Implicación: x <= z
            if req_edges:
                for e in req_edges:
                    solver.Add(x[i, s] <= z[e])

            # Objetivo
            coef_x = int((suitability[s][i] * W_vals[s] * areas[i] * SCALE_SCORE) / SCALE_AREA)
            obj_expr += x[i, s] * coef_x
            
            added = 3.0 - suitability[s][i]
            if added > 0:
                coef_y = int((added * W_vals[s] * areas[i] * SCALE_SCORE) / SCALE_AREA)
                obj_expr += y[i, s] * (coef_y - 1)
            else:
                obj_expr += y[i, s] * -1
                
            active_area_exprs[s] += x[i, s] * areas[i]

    # Costes Z
    for e, var_z in z.items():
        cost_val = edges_unique_cost.get(e, 0)
        cost_expr += var_z * cost_val
        obj_expr += var_z * -1

    # Presupuesto
    solver.Add(cost_expr <= budget_int)

    # Equidad
    total_active = sum(active_area_exprs)
    min_pct = [5, 20, 5, 15]
    for s in range(n_species):
        solver.Add(active_area_exprs[s] * 100 >= min_pct[s] * total_active)

    solver.Maximize(obj_expr)

    # --- EJECUCIÓN ---
    print(f"🧠 Buscando óptimo...")
    status = solver.Solve()
    total_time = time.time() - solver_start_time

    if status in [pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE]:
        obj_val = solver.Objective().Value()
        
        # GAP Calculation
        gap = 0.0
        if status == pywraplp.Solver.OPTIMAL:
            status_str = "OPTIMAL"
        else:
            status_str = "FEASIBLE"
            try:
                best_bound = solver.Objective().BestBound()
                if best_bound != 0:
                    gap = abs(best_bound - obj_val) / abs(best_bound) * 100
            except: gap = -1.0
            
        print("\n" + "="*40)
        print(f"📊 RESULTADO HIGHS")
        print(f"   ⏱️ Tiempo Total:       {total_time:.2f} s")
        print(f"   🎯 Score (Objetivo):   {obj_val / SCALE_SCORE:.2f}")
        print(f"   📉 GAP Final:          {gap:.2f}%")
        print(f"   🚦 Estado:             {status_str}")
        print("="*40)
        
        # Construir CSV
        res = []
        
        # Mapeo rápido de IDs para recuperar nombres
        # (Asumimos que data['ids'] está en orden 0..N-1)
        node_ids = data['ids']
        
        # Pre-calcular qué aristas Z están activas para ir rápido
        # Guardamos en un set las tuplas (u, v) activas
        active_z_edges = set()
        for e, var in z.items():
            # En pywraplp se usa solution_value()
            if var.solution_value() > 0.5:
                active_z_edges.add(e)

        for i in range(n_nodes):
            row = {'grid_id': node_ids[i]}
            
            # 1. Estado de las especies
            for s_idx, s_name in enumerate(S_LIST):
                if paths_requirements[i][s_idx] is not None:
                    val_x = x[i, s_idx].solution_value()
                    val_y = y[i, s_idx].solution_value()
                    row[f'active_{s_name}'] = 1 if val_x > 0.5 else 0
                    row[f'invest_{s_name}'] = 1 if val_y > 0.5 else 0
                else:
                    row[f'active_{s_name}'] = 0
                    row[f'invest_{s_name}'] = 0
            
            # 2. Recuperar Corredores (AHORA SÍ)
            # Miramos los vecinos en la lista de adyacencia
            my_corridors = []
            # data['adj'] es una lista de listas: adj[u] = [(v, cost), ...]
            for (v_idx, _) in data['adj'][i]:
                # La arista en Z siempre está guardada como (min, max)
                edge_key = tuple(sorted((i, v_idx)))
                
                if edge_key in active_z_edges:
                    my_corridors.append(node_ids[v_idx])
            
            row['corridors'] = ";".join(my_corridors)
            res.append(row)

        return pd.DataFrame(res)
    else:
        print(f"❌ Fallo: {status} ({total_time:.2f} s)")
        return None

# ==========================================
# 4. EJECUCIÓN
# ==========================================
if __name__ == "__main__":
    CSV_FILE = 'final_dataset.csv'
    GEOJSON_FILE = 'final_dataset.geojson'
    PRESUPUESTO = 1000
    
    # 1. Pre-procesar
    datos = preprocess_data(CSV_FILE, GEOJSON_FILE, pruning_budget_k=PRESUPUESTO)
    
    # 2. Resolver
    df_final = solve_conservation_highs(datos, budget=PRESUPUESTO, time_limit_sec=600)
    
    if df_final is not None:
        outfile = "solucion_highs_final.csv"
        df_final.to_csv(outfile, index=False)
        print(f"💾 Solución guardada en: {outfile}")