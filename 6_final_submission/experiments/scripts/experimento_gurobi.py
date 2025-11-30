import math
import json
import pandas as pd
import numpy as np
import networkx as nx
import time
import gurobipy as gp
from gurobipy import GRB

# ==========================================
# 1. FUNCIONES AUXILIARES
# ==========================================
# (Idéntico a tu código original)
def haversine_distance(lon1, lat1, lon2, lat2):
    R = 6371.0 
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    return R * c

# ==========================================
# 2. PRE-PROCESAMIENTO
# ==========================================
# (Idéntico a tu código original)
def preprocess_data(csv_path, geojson_path, pruning_budget_k=None):
    print(f"⚙️ [1/2] Pre-procesando datos (Grafo + Dijkstra)...")
    start_time = time.time()
    
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

    cutoff_value = None
    if pruning_budget_k:
        median_cost = np.median(all_edge_costs) if all_edge_costs else 1
        budget_int = int(pruning_budget_k * SCALE_COST)
        limit_pct = budget_int * 0.15 
        safety_floor = median_cost * 6
        cutoff_value = max(limit_pct, safety_floor)

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
# 3. SOLVER (GUROBI)
# ==========================================

def solve_conservation_gurobi(data, budget=500, time_limit_sec=600):
    print(f"\n⚡ [2/2] INICIANDO SOLVER GUROBI (Presupuesto: {budget} k€)")
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

    # Crear modelo Gurobi
    m = gp.Model("Conservation_MIP")
    
    # Parámetros Gurobi (Equivalentes a CP-SAT)
    m.Params.TimeLimit = float(time_limit_sec)
    m.Params.Threads = 8
    m.Params.OutputFlag = 0  # Silenciar log (similar a log_search_progress=False)
    # m.Params.MIPGap = 0.01 # Descomentar si quisieras tolerancia 1%

    # --- VARIABLES ---
    # Gurobi usa addVar/addVars. Usamos binarias (vtype=GRB.BINARY)
    x = {}
    y = {}
    stress = {}
    
    for i in range(n_nodes):
        stress[i] = m.addVar(vtype=GRB.BINARY, name=f"stress_{i}")
        for s in range(n_species):
            x[i, s] = m.addVar(vtype=GRB.BINARY, name=f"x_{i}_{s}")
            y[i, s] = m.addVar(vtype=GRB.BINARY, name=f"y_{i}_{s}")

    # Z: Variables de arista
    used_edges = set()
    for i in range(n_nodes):
        for s in range(n_species):
            reqs = paths_requirements[i][s]
            if reqs: used_edges.update(reqs)
            
    z = {}
    for e in used_edges:
        z[e] = m.addVar(vtype=GRB.BINARY, name=f"z_{e}")

    # --- RESTRICCIONES ---
    # Objetivos parciales
    obj_terms = []
    cost_expr = gp.LinExpr()
    
    active_area_terms = [gp.LinExpr() for _ in range(n_species)]
    PENALTY_STRESS = int(350.0 * SCALE_SCORE)

    for i in range(n_nodes):
        # Stress: x[i,1] + x[i,2] <= 1
        m.addConstr(x[i, 1] + x[i, 2] <= 1)
        # Stress penalty logic: stress >= x_martes + x_conejo - 1
        m.addConstr(stress[i] >= x[i, 1] + x[i, 3] - 1)
        obj_terms.append(stress[i] * (-PENALTY_STRESS))

        for s in range(n_species):
            req_edges = paths_requirements[i][s]
            
            # Poda
            if req_edges is None:
                m.addConstr(x[i, s] == 0)
                continue

            # Inversión: y <= x
            m.addConstr(y[i, s] <= x[i, s])
            
            # Si no es original, x <= y (Inversión obligatoria)
            if not has_spec[s][i]: 
                m.addConstr(x[i, s] <= y[i, s])
            
            # Coste Adaptación
            cost_expr += y[i, s] * cost_adapt[s][i]

            # Implicación de Ruta: x[i,s] -> z[e]  =>  x[i,s] <= z[e]
            if req_edges:
                for e in req_edges:
                    m.addConstr(x[i, s] <= z[e])

            # Objetivo Score
            coef_x = int((suitability[s][i] * W_vals[s] * areas[i] * SCALE_SCORE) / SCALE_AREA)
            obj_terms.append(x[i, s] * coef_x)
            
            # Objetivo Restauración
            added = 3.0 - suitability[s][i]
            if added > 0:
                coef_y = int((added * W_vals[s] * areas[i] * SCALE_SCORE) / SCALE_AREA)
                obj_terms.append(y[i, s] * (coef_y - 1))
            else:
                obj_terms.append(y[i][s] * -1)
                
            active_area_terms[s] += x[i, s] * areas[i]

    # Costes Z
    for e, var_z in z.items():
        cost_val = edges_unique_cost.get(e, 0)
        cost_expr += var_z * cost_val
        obj_terms.append(var_z * -1) # Epsilon negativo

    # Presupuesto
    if budget < 100_000: 
        m.addConstr(cost_expr <= budget_int)

    # Equidad
    # Total Active = sum(Area(s))
    total_active = gp.quicksum(active_area_terms)
    min_pct = [5, 20, 5, 15]
    
    for s in range(n_species):
        # area(s) * 100 >= min_pct * total_active
        m.addConstr(active_area_terms[s] * 100 >= min_pct[s] * total_active)

    # Función Objetivo
    m.setObjective(gp.quicksum(obj_terms), GRB.MAXIMIZE)

    # Ejecución
    print("🧠 Buscando óptimo...")
    m.optimize()
    
    total_time = time.time() - solver_start_time

    # Resultados
    if m.Status in [GRB.OPTIMAL, GRB.TIME_LIMIT, GRB.SUBOPTIMAL]:
        # Verificar si hay solución factible (solcount > 0)
        if m.SolCount > 0:
            obj_val = m.ObjVal
            best_bound = m.ObjBound
            gap = abs(best_bound - obj_val) / abs(best_bound) * 100 if best_bound != 0 else 0.0
                
            # Recuperar estado (Status Code a String)
            status_str = "OPTIMAL" if m.Status == GRB.OPTIMAL else "FEASIBLE/TIME_LIMIT"
            
            # Preparar datos salida
            res = []
            for i in range(n_nodes):
                row = {'grid_id': data['ids'][i]}
                for s_idx, s_name in enumerate(S_LIST):
                    if paths_requirements[i][s_idx] is not None:
                        # Gurobi: .X para obtener valor de variable
                        # (> 0.5 para redondear binario)
                        row[f'active_{s_name}'] = 1 if x[i, s_idx].X > 0.5 else 0
                        row[f'invest_{s_name}'] = 1 if y[i, s_idx].X > 0.5 else 0
                    else:
                        row[f'active_{s_name}'] = 0
                        row[f'invest_{s_name}'] = 0
                res.append(row)
                
            df_sol = pd.DataFrame(res)
            
            metrics = {
                'score': obj_val / SCALE_SCORE,
                'gap': gap,
                'time': total_time,
                'status': status_str
            }
            return df_sol, metrics
        else:
            print(f"❌ Fallo: No integer solution found ({total_time:.2f} s)")
            return None, None
    else:
        print(f"❌ Fallo: Gurobi Status {m.Status}")
        return None, None

# ==========================================
# 4. EJECUCIÓN Y REGISTRO
# ==========================================
if __name__ == "__main__":
    CSV_FILE = 'final_dataset.csv'
    GEOJSON_FILE = 'final_dataset.geojson'
    PRESUPUESTO = 500
    N_RUNS = 30
    
    # 1. Pre-procesar
    datos = preprocess_data(CSV_FILE, GEOJSON_FILE, pruning_budget_k=PRESUPUESTO)
    
    # 2. Bucle de 30 Ejecuciones
    resultados_globales = []
    
    print(f"🚀 Iniciando análisis de robustez con GUROBI ({N_RUNS} ejecuciones)...")
    
    for i in range(N_RUNS):
        print(f"\n🔄 Ejecución {i+1}/{N_RUNS} ...")
        
        # IMPORTANTE: Llamar a la función de Gurobi
        df_sol, metrics = solve_conservation_gurobi(datos, budget=PRESUPUESTO)
        
        if metrics is not None:
            metrics['run_id'] = i
            resultados_globales.append(metrics)
            
            # Opcional: Guardar
            df_res = pd.DataFrame(resultados_globales)
            df_res.to_csv('metricas_robustez_gurobi.csv', index=False)

    # 3. Análisis de Resultados
    df_res = pd.DataFrame(resultados_globales)
    
    print("\n" + "="*40)
    print("📊 RESUMEN ESTADÍSTICO (GUROBI - 30 RUNS)")
    print("="*40)
    print(df_res[['score', 'gap', 'time']].describe().round(2))