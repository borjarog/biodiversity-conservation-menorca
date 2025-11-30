import os
import time
import json
import math
import numpy as np
import pandas as pd
import networkx as nx
from ortools.sat.python import cp_model
from ortools.sat.python.cp_model import CpSolverSolutionCallback
from typing import Optional


class ProgressPrinter(CpSolverSolutionCallback):
    # Visual progress callback for CP-SAT solver

    def __init__(self, start_time):
        CpSolverSolutionCallback.__init__(self)
        self.start_time = start_time
        self.solution_count = 0

    def on_solution_callback(self):
        current_time = time.time()
        elapsed = current_time - self.start_time
        obj = self.ObjectiveValue()
        self.solution_count += 1
        print(
            f"\r [t={elapsed:.1f}s] Solution #{self.solution_count} | Score: {obj:.2f} ...",
            end="",
            flush=True,
        )


class MenorcaSolver:
    def __init__(self, csv_path: str, geojson_path: str):
        self.start_time = time.time()
        self.csv_path = csv_path
        self.geojson_path = geojson_path

        print(
            f" [{self._elapsed()}s] Loading raw data from: {os.path.basename(csv_path)}"
        )
        self.df_raw = pd.read_csv(csv_path)

        # Internal state
        self.df_processed = None
        self.graph_data = None
        self.solution = None

        # Ecological Settings
        self.S_LIST = ["atelerix", "martes", "eliomys", "oryctolagus"]  # Species list
        self.W_VALS = [1.0, 1.0, 2.0, 1.5]  # Species weights

        # Biological Equity [Atelerix, Martes, Eliomys, Oryctolagus]
        self.EQUITY_MIN = [5, 20, 5, 15]  # Minimum percentages
        self.EQUITY_MAX = [30, 60, 30, 50]  # Maximum percentages

        # Scaling
        self.SCALE_COST = 1000
        self.SCALE_SCORE = 10
        self.SCALE_AREA = 100

        self.PENALTY_STRESS = int(350.0 * self.SCALE_SCORE)
        self.RESTORED_Q = 3.0

    def _elapsed(self):
        # Simple elapsed time helper
        return f"{time.time() - self.start_time:.2f}"

    @staticmethod
    def haversine_distance(lon1, lat1, lon2, lat2):
        # Haversine formula to calculate distance between two lat/lon points
        R = 6371.0
        dlat = math.radians(lat2 - lat1)
        dlon = math.radians(lon2 - lon1)
        a = (
            math.sin(dlat / 2) ** 2
            + math.cos(math.radians(lat1))
            * math.cos(math.radians(lat2))
            * math.sin(dlon / 2) ** 2
        )
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        return R * c

    def update_coordinates(self):
        # Updates the dataframe with real-world coordinates from GeoJSON
        print(
            f" [{self._elapsed()}s] Loading real geometry from: {os.path.basename(self.geojson_path)} ..."
        )
        if not os.path.exists(self.geojson_path):
            raise FileNotFoundError(f"Not found: {self.geojson_path}")

        with open(self.geojson_path, "r") as f:
            data = json.load(f)

        centroids = {}
        for feature in data["features"]:
            gid = feature["properties"]["grid_id"]
            geom = feature["geometry"]
            coords = []
            if geom["type"] == "Polygon":
                coords = geom["coordinates"][0]
            elif geom["type"] == "MultiPolygon":
                coords = geom["coordinates"][0][0]

            if coords:
                lon = sum(p[0] for p in coords) / len(coords)
                lat = sum(p[1] for p in coords) / len(coords)
                centroids[gid] = (lon, lat)

        self.df_raw["real_lon"] = self.df_raw["grid_id"].map(
            lambda x: centroids.get(x, (0, 0))[0]
        )
        self.df_raw["real_lat"] = self.df_raw["grid_id"].map(
            lambda x: centroids.get(x, (0, 0))[1]
        )
        print(f" [{self._elapsed()}s] Coordinates updated.")
        self.df_processed = self.df_raw.copy()

    def preprocess_with_pruning(self, pruning_budget_k: float):
        # Preprocesses data and builds graph with pruning based on budget
        t0 = time.time()
        if self.df_processed is None:
            self.update_coordinates()

        df = self.df_processed
        print(
            f" [{self._elapsed()}s] Starting preprocessing with pruning budget: {pruning_budget_k} kEUR..."
        )

        # 1. Basic Data
        ids = df["grid_id"].values
        n_nodes = len(ids)
        id_map = {original_id: i for i, original_id in enumerate(ids)}
        n_species = len(self.S_LIST)

        areas = (df["cell_area_km2"].values * self.SCALE_AREA).astype(int).tolist()
        suitability = [df[f"suitability_{s}"].values.tolist() for s in self.S_LIST]
        cost_adapt = [
            (df[f"cost_adaptation_{s}"].values * self.SCALE_COST).astype(int).tolist()
            for s in self.S_LIST
        ]

        col_map = {
            "atelerix": "has_atelerix_algirus",
            "martes": "has_martes_martes",
            "eliomys": "has_eliomys_quercinus",
            "oryctolagus": "has_oryctolagus_cuniculus",
        }
        has_spec = [df[col_map[s]].values.astype(bool).tolist() for s in self.S_LIST]
        lons = df["real_lon"].values
        lats = df["real_lat"].values

        # 2. Build Physical Graph
        print(f"   [{self._elapsed()}s] Building real-cost graph...")
        G = nx.Graph()
        cost_corr_base = df["cost_corridor"].values
        neighbors_series = df["neighbors"].astype(str).str.split(r"[;,]")

        adj = [[] for _ in range(n_nodes)]
        edges_unique = {}
        all_edge_costs = []

        for u_idx, n_list in enumerate(neighbors_series):
            if not isinstance(n_list, list):
                continue
            for v_str in n_list:
                if v_str in id_map:
                    v_idx = id_map[v_str]
                    if u_idx < v_idx:
                        dist_km = self.haversine_distance(
                            lons[u_idx], lats[u_idx], lons[v_idx], lats[v_idx]
                        )
                        avg_friction = (
                            cost_corr_base[u_idx] + cost_corr_base[v_idx]
                        ) / 2
                        cost = int(avg_friction * dist_km * self.SCALE_COST)

                        G.add_edge(u_idx, v_idx, weight=cost)
                        adj[u_idx].append((v_idx, cost))
                        adj[v_idx].append((u_idx, cost))
                        edges_unique[(u_idx, v_idx)] = cost
                        all_edge_costs.append(cost)

        # 3. Pruning
        median_cost = np.median(all_edge_costs) if all_edge_costs else 1
        budget_int = int(pruning_budget_k * self.SCALE_COST)
        limit_pct = budget_int * 0.15
        safety_floor = median_cost * 6
        cutoff_value = max(limit_pct, safety_floor)
        print(
            f"   [{self._elapsed()}s] Pruning Activated: Ignoring routes > {cutoff_value / self.SCALE_COST:.2f} kEUR"
        )

        # 4. Dijkstra
        print(f"   [{self._elapsed()}s] Calculating all viable routes (Dijkstra)...")
        paths_requirements = [[None for _ in range(n_species)] for _ in range(n_nodes)]

        for s_idx in range(n_species):
            sources = [i for i, val in enumerate(has_spec[s_idx]) if val]
            if not sources:
                continue

            # Pre-compute shortest paths from all source nodes
            dists, paths = nx.multi_source_dijkstra(
                G, sources, weight="weight", cutoff=cutoff_value
            )

            for target, path_nodes in paths.items():
                if target in sources:
                    paths_requirements[target][s_idx] = []
                    continue
                edge_list = []
                for k in range(len(path_nodes) - 1):
                    u, v = sorted((path_nodes[k], path_nodes[k + 1]))
                    edge_list.append((u, v))
                paths_requirements[target][s_idx] = edge_list

        self.graph_data = {
            "ids": ids,
            "n_nodes": n_nodes,
            "id_map": id_map,
            "areas": areas,
            "suitability": suitability,
            "cost_adapt": cost_adapt,
            "has_spec": has_spec,
            "adj": adj,
            "edges_unique": edges_unique,
            "paths_requirements": paths_requirements,
        }
        print(
            f" [{self._elapsed()}s] Preprocessing ready. (Duration: {time.time() - t0:.2f}s)"
        )

    def solve(
        self, budget_k_euro: float, time_limit_sec: int = 600
    ) -> Optional[pd.DataFrame]:
        # Main CP-SAT solver method
        t_solve_start = time.time()

        if self.graph_data is None:
            self.preprocess_with_pruning(budget_k_euro)

        data = self.graph_data
        print(
            f"\n [{self._elapsed()}s] STARTING CP-SAT SOLVER | Budget: {budget_k_euro} kEUR"
        )

        budget_int = int(budget_k_euro * self.SCALE_COST)
        model = cp_model.CpModel()
        n_nodes = data["n_nodes"]
        n_species = len(self.S_LIST)

        # Variables
        # x[i][s]: Cell i active for species s
        x = [
            [model.NewBoolVar(f"x_{i}_{s}") for s in range(n_species)]
            for i in range(n_nodes)
        ]

        # y[i][s]: Cell i adapted for species s
        y = [
            [model.NewBoolVar(f"y_{i}_{s}") for s in range(n_species)]
            for i in range(n_nodes)
        ]

        # stress[i]: Cell i has conflict stress
        stress = [model.NewBoolVar(f"strs_{i}") for i in range(n_nodes)]

        used_edges = set()
        for i in range(n_nodes):
            for s in range(n_species):
                reqs = data["paths_requirements"][i][s]
                if reqs:
                    used_edges.update(reqs)

        # z[e]: Edge e used in any corridor
        z = {}
        for e in used_edges:
            z[e] = model.NewBoolVar(f"z_{e}")

        # Constraints
        obj_terms = []  # Objective terms
        cost_terms = []
        active_area_terms = [[] for _ in range(n_species)]

        # Inter-species conflict and objective terms
        for i in range(n_nodes):
            model.Add(x[i][1] + x[i][2] <= 1)
            model.Add(stress[i] >= x[i][1] + x[i][3] - 1)
            obj_terms.append(stress[i] * (-self.PENALTY_STRESS))

            for s in range(n_species):
                req_edges = data["paths_requirements"][i][s]
                if req_edges is None:
                    model.Add(x[i][s] == 0)
                    continue

                # If cell does not have the species, cannot be active without adaptation
                model.Add(y[i][s] <= x[i][s])
                if not data["has_spec"][s][i]:
                    model.Add(x[i][s] <= y[i][s])
                cost_terms.append(y[i][s] * data["cost_adapt"][s][i])

                # If active, must build corridors
                if req_edges:
                    for e in req_edges:
                        model.AddImplication(x[i][s], z[e])

                coef_x = int(
                    (
                        data["suitability"][s][i]
                        * self.W_VALS[s]
                        * data["areas"][i]
                        * self.SCALE_SCORE
                    )
                    / self.SCALE_AREA
                )
                obj_terms.append(x[i][s] * coef_x)

                added = self.RESTORED_Q - data["suitability"][s][i]
                if added > 0:
                    coef_y = int(
                        (added * self.W_VALS[s] * data["areas"][i] * self.SCALE_SCORE)
                        / self.SCALE_AREA
                    )
                    obj_terms.append(y[i][s] * (coef_y - 1))
                else:
                    obj_terms.append(y[i][s] * -1)

                active_area_terms[s].append(x[i][s] * data["areas"][i])

        # Corridor costs and objective terms
        for e, z_var in z.items():
            cost_terms.append(z_var * data["edges_unique"][e])
            obj_terms.append(z_var * -1)

        # Budget Constraint
        if budget_k_euro < 100_000:
            model.Add(sum(cost_terms) <= budget_int)

        # Equity Constraints
        total_active = sum(sum(active_area_terms[s]) for s in range(n_species))
        for s in range(n_species):
            s_sum = sum(active_area_terms[s])
            model.Add(s_sum * 100 >= self.EQUITY_MIN[s] * total_active)
            model.Add(s_sum * 100 <= self.EQUITY_MAX[s] * total_active)

        model.Maximize(sum(obj_terms))

        # Solver Config
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = time_limit_sec
        solver.parameters.num_search_workers = 12
        solver.parameters.log_search_progress = False

        print(f"[{self._elapsed()}s] Sending to solver (Limit: {time_limit_sec}s)...")
        print("   (Updating progress in real-time...)")

        # Callback for loading bar
        callback = ProgressPrinter(time.time())
        status = solver.Solve(model, callback)
        print("\n")

        t_total_solve = time.time() - t_solve_start

        # Results
        status_name = solver.StatusName(status)
        obj_val = solver.ObjectiveValue() / self.SCALE_SCORE
        best_bound = solver.BestObjectiveBound() / self.SCALE_SCORE

        gap = 0.0
        if best_bound > 0:
            gap = abs(best_bound - obj_val) / abs(best_bound) * 100

        print("\n" + "=" * 40)
        print(f" TOTAL TIME: {t_total_solve:.2f} seconds")
        if status == cp_model.OPTIMAL:
            print(f" OPTIMAL SOLUTION PROVEN")
        elif status == cp_model.FEASIBLE:
            print(f" FEASIBLE SOLUTION (Gap: {gap:.2f}%)")
        else:
            print(f" FAILURE: {status_name}")
            return None

        print(f" Final Score: {obj_val:.2f}")
        print("=" * 40 + "\n")

        self.solution = self._format_solution(solver, x, y, z, data)
        return self.solution

    def _format_solution(self, solver, x, y, z, data):
        res = []
        ids = data["ids"]
        for i in range(data["n_nodes"]):
            row = {"grid_id": ids[i]}
            for s_idx, s_name in enumerate(self.S_LIST):
                if data["paths_requirements"][i][s_idx] is not None:
                    row[f"active_{s_name}"] = solver.Value(x[i][s_idx])
                    row[f"invest_{s_name}"] = solver.Value(y[i][s_idx])
                else:
                    row[f"active_{s_name}"] = 0
                    row[f"invest_{s_name}"] = 0

            cors = []
            for v_idx, _ in data["adj"][i]:
                u, v = sorted((i, v_idx))
                edge = (u, v)
                if edge in z and solver.Value(z[edge]):
                    cors.append(ids[v_idx])
            row["corridors"] = ";".join(cors)
            res.append(row)
        return pd.DataFrame(res)

    def audit_results(self):
        # Audits and prints cost breakdown of the solution
        if self.solution is None:
            print(" Warning: No solution loaded.")
            return

        print(f" COST AUDIT")
        print("-" * 50)

        df_sol = self.solution
        df_base = self.df_processed
        df_merged = pd.merge(
            df_sol, df_base, on="grid_id", how="left", suffixes=("", "_base")
        )

        # 1. Adaptation
        total_adapt = 0
        for s in self.S_LIST:
            inv_col = f"invest_{s}"
            cost_col = f"cost_adaptation_{s}"
            if inv_col in df_merged.columns:
                sub = df_merged[df_merged[inv_col] == 1]
                cost = sub[cost_col].sum()
                total_adapt += cost
                print(
                    f"  > {s.capitalize():<12}: {len(sub):4d} cells | {cost:8.2f} kEUR"
                )

        # 2. Corridors
        print("-" * 50)
        coords = df_base.set_index("grid_id")[["real_lon", "real_lat"]].to_dict("index")
        corr_costs = df_base.set_index("grid_id")["cost_corridor"].to_dict()

        total_corr = 0
        total_km = 0
        processed = set()

        for idx, row in df_sol.iterrows():
            u = row["grid_id"]
            neighbors = str(row.get("corridors", "")).replace(",", ";").split(";")
            for v in neighbors:
                v = v.strip()
                if v and v in coords:
                    edge = tuple(sorted((u, v)))
                    if edge not in processed:
                        processed.add(edge)
                        d = self.haversine_distance(
                            coords[u]["real_lon"],
                            coords[u]["real_lat"],
                            coords[v]["real_lon"],
                            coords[v]["real_lat"],
                        )
                        friction = (corr_costs[u] + corr_costs[v]) / 2
                        total_corr += friction * d
                        total_km += d

        print(f"  > Corridors:     {total_km:8.2f} km")
        print(f"  > Infra Cost:    {total_corr:8.2f} kEUR")
        print("-" * 50)
        print(f" FINAL TOTAL:     {total_adapt + total_corr:10.2f} kEUR")
