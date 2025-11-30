import math
import json
import pandas as pd
import numpy as np
import networkx as nx
import time
from ortools.linear_solver import pywraplp
import multiprocessing
import os # Necesario para manejar procesos

# ==========================================
# 1. FUNCIONES AUXILIARES
# ==========================================
def haversine_distance(lon1, lat1, lon2, lat2):
    """Calcula la distancia Haversine entre dos puntos (lon, lat) en km."""
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
    """Carga y pre-procesa todos los datos para el MIP."""
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
        # Desactivamos el cutoff para multi_source_dijkstra ya que no está implementado
        # en la versión de networkx que usamos
        try:
             # Usamos la versión sin cutoff que es compatible con la firma
            dists, paths = nx.multi_source_dijkstra(G, sources, weight="weight") 
        except TypeError:
            # En caso de necesitar un cutoff, usaríamos una implementación a medida
            dists, paths = nx.multi_source_dijkstra(G, sources, weight="weight") 
            # Filtrado manual si es necesario:
            # paths = {k: v for k, v in paths.items() if dists.get(k, float('inf')) <= cutoff_value}
        
        for tgt,path_nodes in paths.items():
            if tgt in sources:
                paths_requirements[tgt][s_idx] = []
                continue
            # Verifica la distancia para aplicar el cutoff si la versión de NX no lo hizo
            if cutoff_value is not None and dists.get(tgt, float('inf')) > cutoff_value:
                continue

            edge_list = []
            for k in range(len(path_nodes)-1):
                u,v = sorted((path_nodes[k], path_nodes[k+1]))
                edge_list.append((u,v))
            paths_requirements[tgt][s_idx] = edge_list

    print(f"✅ Pre-procesamiento en {time.time()-start_time:.2f}s")
    # Es crucial que data sea 'picklable' para multiprocessing, un diccionario lo es.
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
    # El parámetro 'seed' se mantiene en la firma pero ya no se utiliza
    print(f"\n⚡ Ejecutando {solver_name} (Run={seed}, PID={os.getpid()})")

    # 1. CREACIÓN DEL SOLVER
    if solver_name == "HIGHS":
        solver = pywraplp.Solver.CreateSolver("HIGHS")
    elif solver_name == "CBC":
        solver = pywraplp.Solver.CreateSolver("CBC")
    elif solver_name == "SCIP":
        solver = pywraplp.Solver.CreateSolver("SCIP")
    else:
        solver = pywraplp.Solver.CreateSolver("SCIP") # Fallback

    if not solver:
        return None, None

    # 2. CONFIGURACIÓN GENERAL (SIN SEMILLA ESPECÍFICA)
    solver.SetNumThreads(1)
    solver.SetTimeLimit(time_limit_sec * 1000)
    # Se elimina la lógica de SetSolverSpecificParametersAsString para evitar errores.
    # El solver usará su semilla aleatoria por defecto (basada en el reloj del sistema).
    
    # DESCOMPACTACIÓN
    n = data['n_nodes']
    ns = len(data['S_LIST'])
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
    x = {} # Conservación de célula i para especie s
    y = {} # Adaptación/mejora de célula i para especie s
    stress = {} # Tensión entre 'martes' y 'oryctolagus'

    used_edges = set()

    for i in range(n):
        stress[i] = solver.IntVar(0,1,"stress_%d"%i)
        for s in range(ns):
            x[i,s] = solver.IntVar(0,1,f"x_{i}_{s}")
            y[i,s] = solver.IntVar(0,1,f"y_{i}_{s}")
            if paths[i][s]:
                used_edges.update(paths[i][s])

    z = {} # Creación de corredor (edge) e
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
                coef_y = int((added * W[s] * areas[i] * SCALE_SCORE) / SCALE_AREA)
                
                obj += (coef_y - 1) * y[i,s] 
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

    # Iniciar el tiempo de resolución
    start_solve_time = time.time()
    status = solver.Solve()
    elapsed = time.time() - start_solve_time

    if status not in (pywraplp.Solver.OPTIMAL, pywraplp.Solver.FEASIBLE):
        print(f"❌ {solver_name} (Run {seed}): Status no óptimo/factible: {status}")
        return None, None
    
    # ... (Lógica de métricas) ...

    if status == pywraplp.Solver.OPTIMAL:
        obj_val = solver.Objective().Value() / SCALE_SCORE
        gap = 0.0
    else: # FEASIBLE o TIMEOUT
        obj_val = solver.Objective().Value() / SCALE_SCORE 
        best_bound = solver.Objective().BestBound()
        if best_bound > 0 and best_bound >= obj_val * SCALE_SCORE:
            gap = (best_bound - obj_val * SCALE_SCORE) / best_bound 
        else:
            gap = 1.0 
            
    if elapsed > time_limit_sec:
        elapsed = time_limit_sec
        status_str = "T/O"
    elif status == pywraplp.Solver.OPTIMAL:
        status_str = "OPT"
    else:
        status_str = "FEAS"

    print(f"✅ {solver_name} (Run {seed}): Score={obj_val:.2f}, Time={elapsed:.2f}s, Status={status_str}")

    return {}, {
        "score": obj_val,
        "time": elapsed,
        "status": status_str,
        "gap": gap, 
        "solver_internal_gap": solver.Objective().BestBound(),
    }


# ==========================================
# 4. FUNCIÓN DEL EXPERIMENTO PARA UN SOLO SOLVER
# ==========================================
def run_solver_experiment(solver_name, data, budget, n_runs, max_time_per_solver):
    """Ejecuta el experimento completo para un único solver."""
    print(f"\n==============================")
    print(f"   🚀 EXPERIMENTO {solver_name}")
    print(f"     (Límite: {max_time_per_solver/60:.0f} min)")
    print("==============================\n")

    results = []
    start_solver_time = time.time()

    # Cada run tiene un tiempo límite de 600 segundos (10 minutos)
    run_time_limit_sec = 600 

    for run in range(n_runs):
        # Verifica el tiempo máximo TOTAL para este solver
        current_elapsed = time.time() - start_solver_time
        if current_elapsed > max_time_per_solver:
            print(f"⏳ Tiempo máximo ({max_time_per_solver/60:.0f} min) alcanzado para {solver_name}. Se detiene después de {run} runs.")
            break

        print(f"🔄 {solver_name}: Run {run+1}/{n_runs}")
        
        # El tiempo límite para la resolución individual se mantiene en 600s
        df_sol, metrics = solve_conservation_mip(
            data,
            budget=budget,
            time_limit_sec=run_time_limit_sec, 
            seed=run, # Mantenemos 'run' solo como ID para el output y el CSV.
            solver_name=solver_name
        )

        if metrics:
            metrics["run_id"] = run
            metrics["solver"] = solver_name
            results.append(metrics)
        else:
             print(f"⚠️ {solver_name}: Run {run+1} falló o no encontró solución.")


    df_res = pd.DataFrame(results)
    # Asegura que el nombre del archivo sea único y descriptivo
    output_filename = f"metricas_{solver_name}_{budget}_{run_time_limit_sec//60}min_multiprocess_noseed.csv"
    df_res.to_csv(output_filename, index=False)

    print(f"\n📊 Resumen de {solver_name}:")
    print(df_res[["score","time"]].describe().round(2))
    print(f"💾 Resultados guardados en {output_filename}")


# ==========================================
# 5. EXPERIMENTO PRINCIPAL CON MULTIPROCESSING
# ==========================================
if __name__ == "__main__":
    # La validación de la compatibilidad con multiprocessing es esencial
    # en Windows o si el script se ejecuta de forma diferente a la habitual.
    multiprocessing.freeze_support()
    
    CSV_FILE = "final_dataset.csv"
    GEOJSON_FILE = "final_dataset.geojson"
    PRESUPUESTO = 500
    N_RUNS = 30
    MAX_SOLVER_TIME = 4500   # 1h 15m (Límite TOTAL por solver)
    
    # 1. PRE-PROCESAMIENTO (Compartido)
    # Se ejecuta una sola vez al inicio.
    data = preprocess_data(CSV_FILE, GEOJSON_FILE, pruning_budget_k=PRESUPUESTO)
    
    # 2. CONFIGURACIÓN DE LOS SOLVERS
    SOLVERS = ["CBC", "SCIP", "HIGHS"]
    processes = []
    
    print("\n\n#####################################################")
    print(f"🚀 INICIANDO EXPERIMENTO PARALELO con {len(SOLVERS)} Solvers.")
    print(f"⏰ Límite de tiempo global por solver: {MAX_SOLVER_TIME}s")
    print("#####################################################\n")
    
    global_start_time = time.time()

    # 3. INICIO DE PROCESOS PARALELOS
    for solver_name in SOLVERS:
        # Crea un nuevo proceso para cada solver
        process_args = (solver_name, data, PRESUPUESTO, N_RUNS, MAX_SOLVER_TIME)
        p = multiprocessing.Process(target=run_solver_experiment, args=process_args)
        processes.append(p)
        p.start()
        print(f"✅ Proceso para {solver_name} iniciado con PID: {p.pid}")

    # 4. ESPERAR A QUE TODOS LOS PROCESOS TERMINEN
    
    for p in processes:
        p.join() 
    
    total_elapsed = time.time() - global_start_time
    
    print("\n#####################################################")
    print("🏁 EXPERIMENTO PARALELO FINALIZADO")
    print(f"⏳ Tiempo TOTAL de ejecución: {total_elapsed:.2f} segundos.")
    print("#####################################################")

    print("\nSe han generado los archivos CSV con los resultados de cada solver.")