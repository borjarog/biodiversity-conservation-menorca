import math
import json
import pandas as pd
import numpy as np
import networkx as nx
import time
from ortools.linear_solver import pywraplp
import os

# ==========================================
# 1. FUNCIONES AUXILIARES
# ==========================================
def haversine_distance(lon1, lat1, lon2, lat2):
    R = 6371.0
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1))*math.cos(math.radians(lat2))*math.sin(dlon/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return R * c

# ==========================================
# 2. PRE-PROCESAMIENTO
# ==========================================
def preprocess_data(csv_path, geojson_path, pruning_budget_k=None):
    print("⚙️ Pre-procesando datos...")
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
    
    areas = (df['cell_area_km2'] * SCALE_AREA).astype(int).tolist()
    suitability = [df[f'suitability_{s}'].tolist() for s in S_LIST]
    cost_adapt = [(df[f'cost_adaptation_{s}'] * SCALE_COST).astype(int).tolist() for s in S_LIST]

    col_map = {
        'atelerix':'has_atelerix_algirus',
        'martes':'has_martes_martes',
        'eliomys':'has_eliomys_quercinus',
        'oryctolagus':'has_oryctolagus_cuniculus'
    }
    has_spec = [df[col_map[s]].astype(bool).tolist() for s in S_LIST]

    with open(geojson_path, "r") as f:
        geo = json.load(f)

    centroids = {}
    for feat in geo['features']:
        gid = feat['properties']['grid_id']
        geom = feat['geometry']
        if geom['type'] == 'Polygon':
            coords = geom['coordinates'][0]
        else:
            coords = geom['coordinates'][0][0]
        lon = sum(p[0] for p in coords)/len(coords)
        lat = sum(p[1] for p in coords)/len(coords)
        centroids[gid] = (lon, lat)
    
    lons = np.array([centroids[i][0] for i in ids])
    lats = np.array([centroids[i][1] for i in ids])

    G = nx.Graph()
    cost_corr = df['cost_corridor'].values
    neighbors_series = df['neighbors'].astype(str).str.split(r"[;,]")

    edges_unique_cost = {}
    all_edge_costs = []

    for u_idx, lst in enumerate(neighbors_series):
        if not isinstance(lst, list):
            continue
        for v_str in lst:
            if v_str in id_map:
                v_idx = id_map[v_str]
                if u_idx < v_idx:
                    dist_km = haversine_distance(lons[u_idx], lats[u_idx], lons[v_idx], lats[v_idx])
                    avg_frict = (cost_corr[u_idx] + cost_corr[v_idx]) / 2
                    cost = int(avg_frict * dist_km * SCALE_COST)
                    G.add_edge(u_idx, v_idx, weight=cost)
                    edges_unique_cost[(u_idx,v_idx)] = cost
                    all_edge_costs.append(cost)

    cutoff_value = None
    if pruning_budget_k:
        median_cost = np.median(all_edge_costs) if all_edge_costs else 1
        budget_int = pruning_budget_k * SCALE_COST
        cutoff_value = max(budget_int * 0.15, median_cost * 6)

    paths_requirements = [[None]*n_species for _ in range(n_nodes)]

    for s_idx in range(n_species):
        sources = [i for i,b in enumerate(has_spec[s_idx]) if b]
        if not sources:
            continue
        dists, paths = nx.multi_source_dijkstra(G, sources, weight="weight")
        for tgt,path_nodes in paths.items():
            if tgt in sources:
                paths_requirements[tgt][s_idx] = []
                continue
            if cutoff_value is not None and dists.get(tgt,float('inf'))>cutoff_value:
                continue
            edge_list = []
            for k in range(len(path_nodes)-1):
                u,v = sorted((path_nodes[k], path_nodes[k+1]))
                edge_list.append((u,v))
            paths_requirements[tgt][s_idx] = edge_list

    print(f"✅ Pre-procesamiento en {time.time()-start_time:.2f}s")
    return {
        'ids': ids, 'n_nodes': n_nodes, 'S_LIST': S_LIST, 'W_vals': W_vals,
        'areas': areas, 'suitability': suitability, 'cost_adapt': cost_adapt,
        'has_spec': has_spec, 'edges_unique_cost': edges_unique_cost,
        'paths_requirements': paths_requirements,
        'SCALE_COST': SCALE_COST, 'SCALE_SCORE': SCALE_SCORE, 'SCALE_AREA': SCALE_AREA
    }

# ==========================================
# 3. SOLVER HIGHS 8 THREADS
# ==========================================
def solve_highs(data, budget=1000, time_limit_sec=600):
    print("\n⚡ Ejecutando HIGHS con 8 threads")
    solver = pywraplp.Solver.CreateSolver("HIGHS")
    if not solver:
        print("❌ No se pudo crear el solver HIGHS")
        return None, None

    solver.SetNumThreads(8)          # 8 threads
    solver.SetTimeLimit(time_limit_sec * 1000)  # 10 minutos

    # Desempaquetado simplificado
    n = data['n_nodes']
    ns = len(data['S_LIST'])
    areas = data['areas']
    W = data['W_vals']
    suitability = data['suitability']
    cost_adapt = data['cost_adapt']
    has_spec = data['has_spec']
    edges_cost = data['edges_unique_cost']
    paths = data['paths_requirements']
    SCALE_COST = data['SCALE_COST']
    SCALE_SCORE = data['SCALE_SCORE']
    SCALE_AREA = data['SCALE_AREA']

    budget_int = budget * SCALE_COST

    # VARIABLES
    x = {}
    y = {}
    stress = {}
    used_edges = set()
    for i in range(n):
        stress[i] = solver.IntVar(0,1,f"stress_{i}")
        for s in range(ns):
            x[i,s] = solver.IntVar(0,1,f"x_{i}_{s}")
            y[i,s] = solver.IntVar(0,1,f"y_{i}_{s}")
            if paths[i][s]:
                used_edges.update(paths[i][s])
    z = {e: solver.IntVar(0,1,f"z_{e}") for e in used_edges}

    # OBJETO Y RESTRICCIONES simplificado (igual que antes)
    obj = 0
    cost_expr = 0
    species_area = [0]*ns
    PENALTY_STRESS = int(350*SCALE_SCORE)
    for i in range(n):
        solver.Add(x[i,1]+x[i,2] <= 1)
        solver.Add(stress[i] >= x[i,1]+x[i,3]-1)
        obj += -PENALTY_STRESS*stress[i]

        for s in range(ns):
            req = paths[i][s]
            solver.Add(y[i,s]<=x[i,s])
            if not has_spec[s][i]:
                solver.Add(x[i,s]<=y[i,s])
            cost_expr += cost_adapt[s][i]*y[i,s]
            if req:
                for e in req:
                    solver.Add(x[i,s]<=z[e])
            coef_x = int((suitability[s][i]*W[s]*areas[i]*SCALE_SCORE)/SCALE_AREA)
            obj += coef_x*x[i,s]
            added = 3 - suitability[s][i]
            if added>0:
                coef_y = int((added*W[s]*areas[i]*SCALE_SCORE)/SCALE_AREA)
                obj += (coef_y-1)*y[i,s]
            else:
                obj += -1*y[i,s]
            species_area[s] += areas[i]*x[i,s]

    for e,var in z.items():
        obj += -1*var
        cost_expr += edges_cost[e]*var

    solver.Add(cost_expr<=budget_int)
    total_active = sum(species_area)
    min_pct = [5,20,5,15]
    for s in range(ns):
        solver.Add(species_area[s]*100 >= min_pct[s]*total_active)
    solver.Maximize(obj)

    # RESOLVER
    start_time = time.time()
    status = solver.Solve()
    elapsed = time.time()-start_time

    if status not in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE):
        print("❌ HIGHS no encontró solución")
        return None, None

    if status == pywraplp.Solver.OPTIMAL:
        status_str="OPT"
        obj_val = solver.Objective().Value()/SCALE_SCORE
        gap = 0.0
    else:
        status_str="FEAS"
        obj_val = solver.Objective().Value()/SCALE_SCORE
        try:
            best_bound = solver.Objective().BestBound()
            gap = (best_bound - obj_val*SCALE_SCORE)/best_bound
        except:
            gap=1.0

    print(f"✅ HIGHS: Score={obj_val:.2f}, Gap={gap:.2f}, Tiempo={elapsed:.2f}s, Status={status_str}")
    return {}, {"score":obj_val, "gap":gap, "time":elapsed, "status":status_str}


# ==========================================
# 4. EJECUCIÓN SIMPLE
# ==========================================
if __name__ == "__main__":
    CSV_FILE = "final_dataset.csv"
    GEOJSON_FILE = "final_dataset.geojson"
    PRESUPUESTO = 500
    N_RUNS = 30
    TIME_LIMIT = 600  # 10 min por run

    data = preprocess_data(CSV_FILE, GEOJSON_FILE, pruning_budget_k=PRESUPUESTO)

    resultados = []
    for i in range(N_RUNS):
        print(f"\n🔄 Ejecución {i+1}/{N_RUNS}")
        _, metrics = solve_highs(data, budget=PRESUPUESTO, time_limit_sec=TIME_LIMIT)
        if metrics:
            metrics["run_id"]=i
            resultados.append(metrics)

    df_res = pd.DataFrame(resultados)
    df_res.to_csv(f"metricas_highs_8threads_{PRESUPUESTO}.csv", index=False)
    print("\n📊 Resumen final:")
    print(df_res.describe().round(2))
