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
    R = 6371.0 
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1))*math.cos(math.radians(lat2))*math.sin(dlon/2)**2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    return R * c


# ==========================================
# 2. PRE-PROCESAMIENTO (idéntico al tuyo)
# ==========================================
def preprocess_data(csv_path, geojson_path, pruning_budget_k=None):
    print(f"⚙️ [1/2] Pre-procesando datos...")
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
        dists, paths = nx.multi_source_dijkstra(G, sources, weight="weight", cutoff=cutoff_value)
        for tgt,path_nodes in paths.items():
            if tgt in sources:
                paths_requirements[tgt][s_idx] = []
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
# 3. SOLVER UNIFICADO (1 THREAD, 10min)
# ==========================================
def solve_conservation_mip(data, budget=1000, time_limit_sec=600, seed=0, solver_name='CBC'):
    print(f"\n⚡ Ejecutando {solver_name} (Seed={seed})")

    # Crear solver con 1 thread y 10 min por ejecución
    if solver_name == "HIGHS":
        solver = pywraplp.Solver.CreateSolver("HIGHS")
    elif solver_name == "CBC":
        solver = pywraplp.Solver.CreateSolver("CBC")
    else:
        solver = pywraplp.Solver.CreateSolver("SCIP")

    if not solver:
        return None, None

    solver.SetNumThreads(1)
    solver.SetTimeLimit(time_limit_sec * 1000)

    # DESCOMPACTACIÓN
    n = data['n_nodes']
    S = data['S_LIST']
    ns = len(S)
    W = data['W_vals']
    suitability = data['suitability']
    cost_adapt = data['cost_adapt']
    has_spec = data['has_spec']
    areas = data['areas']
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
        stress[i] = solver.IntVar(0,1,"stress_%d"%i)
        for s in range(ns):
            x[i,s] = solver.IntVar(0,1,f"x_{i}_{s}")
            y[i,s] = solver.IntVar(0,1,f"y_{i}_{s}")
            if paths[i][s]:
                used_edges.update(paths[i][s])

    z = {}
    for e in used_edges:
        z[e] = solver.IntVar(0,1,"z_%s_%s"%e)

    # OBJETO Y RESTRICCIONES
    obj = 0
    cost_expr = 0
    species_area = [0]*ns
    PENALTY_STRESS = int(350*SCALE_SCORE)

    for i in range(n):
        solver.Add(x[i,1] + x[i,2] <= 1)
        solver.Add(stress[i] >= x[i,1] + x[i,3] - 1)
        obj += -PENALTY_STRESS * stress[i]

        for s in range(ns):
            req = paths[i][s]

            solver.Add(y[i,s] <= x[i,s])
            if not has_spec[s][i]:
                solver.Add(x[i,s] <= y[i,s])

            cost_expr += cost_adapt[s][i] * y[i,s]

            if req:
                for e in req:
                    solver.Add(x[i,s] <= z[e])

            coef_x = int((suitability[s][i] * W[s] * areas[i] * SCALE_SCORE) / SCALE_AREA)
            obj += coef_x * x[i,s]

            added = 3 - suitability[s][i]
            if added > 0:
                obj += int((added * W[s] * areas[i] * SCALE_SCORE) / SCALE_AREA) * y[i,s]
            else:
                obj += -1 * y[i,s]

            species_area[s] += areas[i] * x[i,s]

    for e,var in z.items():
        obj += -1 * var
        cost_expr += edges_cost[e] * var

    solver.Add(cost_expr <= budget_int)

    total_active = sum(species_area)
    min_pct = [5,20,5,15]
    for s in range(ns):
        solver.Add(species_area[s] * 100 >= min_pct[s] * total_active)

    solver.Maximize(obj)

    status = solver.Solve()
    elapsed = solver.wall_time()/1000

    if status not in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE):
        return None, None

    obj_val = solver.Objective().Value() / SCALE_SCORE

    return {}, {
        "score": obj_val,
        "time": elapsed,
        "status": "OPT" if status==pywraplp.Solver.OPTIMAL else "FEAS",
        "gap": solver.Objective().BestBound(), 
    }


# ==========================================
# 4. EXPERIMENTO (MÁX 1h 15m por SOLVER)
# ==========================================
if __name__ == "__main__":
    CSV_FILE = "final_dataset.csv"
    GEOJSON_FILE = "final_dataset.geojson"
    PRESUPUESTO = 1000
    N_RUNS = 30
    MAX_SOLVER_TIME = 4500   # 1h 15m

    data = preprocess_data(CSV_FILE, GEOJSON_FILE, pruning_budget_k=PRESUPUESTO)
    # lista = ["CBC", "SCIP", "HIGHS"]
    for solver_name in ["SCIP"]:
        print(f"\n==============================")
        print(f"   🚀 EXPERIMENTO {solver_name}")
        print("==============================\n")

        results = []
        start_solver_time = time.time()

        for run in range(N_RUNS):
            if time.time() - start_solver_time > MAX_SOLVER_TIME:
                print(f"⏳ Tiempo máximo alcanzado para {solver_name}. Se detiene.")
                break

            print(f"🔄 Run {run+1}/{N_RUNS}")
            df_sol, metrics = solve_conservation_mip(
                data,
                budget=PRESUPUESTO,
                time_limit_sec=600,
                seed=run,
                solver_name=solver_name
            )

            if metrics:
                metrics["run_id"] = run
                metrics["solver"] = solver_name
                results.append(metrics)

        df_res = pd.DataFrame(results)
        df_res.to_csv(f"metricas_{solver_name}_500_1_thread.csv", index=False)

        print("\n📊 Resumen:")
        print(df_res[["score","time"]].describe().round(2))
