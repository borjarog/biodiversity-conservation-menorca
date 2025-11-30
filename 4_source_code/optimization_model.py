"""
Menorca Biodiversity Conservation Solver (Model V.3)
====================================================

This module implements the definitive optimization pipeline for the Menorca conservation project.
It solves a Mixed-Integer Programming (MIP) problem to select optimal habitat cells and corridors
under budget constraints, ensuring connectivity and species equity.

Key Features:
    - Hybrid 3-stage pipeline: Preprocessing (Graph) -> Optimization (CP-SAT) -> Auditing.
    - Real Geography: Uses Haversine distances and friction costs from a GeoJSON grid.
    - Smart Pruning: Uses multi-source Dijkstra to pre-calculate valid paths and prune the search space.
    - Constraints: Budget, Connectivity (Path Implication), Conflict (Martes/Eliomys), Equity (Area Shares).

Usage:
    Import the class and run the solver:

    solver = MenorcaSolver("data.csv", "map.geojson")
    df_solution = solver.solve(budget_k_euro=1000)
    solver.audit_results()
"""

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
    # Callback to print the progress of the CP-SAT solver to the console in real-time.

    def __init__(self, start_time):
        CpSolverSolutionCallback.__init__(self)
        self.start_time = start_time
        self.solution_count = 0

    def on_solution_callback(self):
        """Called by the solver each time a new feasible solution is found."""
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
    """Main solver class encapsulating the data loading, preprocessing, and optimization logic."""

    def __init__(self, csv_path: str, geojson_path: str):
        """
        Initializes the solver with paths to the input data.

        Args:
            csv_path (str): Path to the CSV file containing grid attributes (suitability, cost, etc.).
            geojson_path (str): Path to the GeoJSON file containing the grid geometry (polygons).
        """
        self.start_time = time.time()
        self.csv_path = csv_path
        self.geojson_path = geojson_path

        print(
            f" [{self._elapsed()}s] Loading raw data from: {os.path.basename(csv_path)}"
        )
        self.df_raw = pd.read_csv(csv_path)

        # Internal state variables
        self.df_processed = None
        self.graph_data = None  # Dictionary to store pre-calculated graph & paths
        self.solution = None  # DataFrame storing the final result

        # --- Ecological Parameters ---
        self.S_LIST = ["atelerix", "martes", "eliomys", "oryctolagus"]  # Species list
        self.W_VALS = [1.0, 1.0, 2.0, 1.5]  # Ecological weights per species

        # Equity Constraints (Min/Max % of total area)
        # Order matches S_LIST: [Atelerix, Martes, Eliomys, Oryctolagus]
        self.EQUITY_MIN = [5, 20, 5, 15]
        self.EQUITY_MAX = [30, 60, 30, 50]

        # --- Solver Scaling Factors (CP-SAT requires integers) ---
        self.SCALE_COST = 1000  # 1 k€ becomes 1000 units
        self.SCALE_SCORE = 10  # Suitability 3.0 becomes 30 units
        self.SCALE_AREA = 100  # Area 0.25 km2 becomes 25 units

        self.PENALTY_STRESS = int(350.0 * self.SCALE_SCORE)
        self.RESTORED_Q = 3.0

    def _elapsed(self):
        """Returns formatted elapsed time string since initialization."""
        return f"{time.time() - self.start_time:.2f}"

    @staticmethod
    def haversine_distance(lon1, lat1, lon2, lat2):
        """
        Calculates the great-circle distance between two points on Earth.

        Args:
            lon1, lat1: Coordinates of point A (in degrees).
            lon2, lat2: Coordinates of point B (in degrees).

        Returns:
            float: Distance in kilometers.
        """
        R = 6371.0  # Earth radius in km
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
        """
        Reads the GeoJSON file to extract the real centroid (lon/lat) for each grid cell.
        Updates self.df_raw with 'real_lon' and 'real_lat' columns.
        """
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

            # Handle Polygon and MultiPolygon geometries
            if geom["type"] == "Polygon":
                coords = geom["coordinates"][0]  # Outer ring
            elif geom["type"] == "MultiPolygon":
                coords = geom["coordinates"][0][0]  # First polygon outer ring

            if coords:
                # Calculate simple centroid average
                lon = sum(p[0] for p in coords) / len(coords)
                lat = sum(p[1] for p in coords) / len(coords)
                centroids[gid] = (lon, lat)

        # Map centroids to DataFrame
        self.df_raw["real_lon"] = self.df_raw["grid_id"].map(
            lambda x: centroids.get(x, (0, 0))[0]
        )
        self.df_raw["real_lat"] = self.df_raw["grid_id"].map(
            lambda x: centroids.get(x, (0, 0))[1]
        )
        print(f" [{self._elapsed()}s] Coordinates updated.")
        self.df_processed = self.df_raw.copy()

    def preprocess_with_pruning(self, pruning_budget_k: float):
        """
        Builds the connectivity graph and pre-calculates shortest paths for all species.
        Applies pruning to discard routes that are too expensive relative to the budget.

        Args: pruning_budget_k (float): The reference budget (k€) used to calculate the cutoff threshold.
        """
        t0 = time.time()
        if self.df_processed is None:
            self.update_coordinates()

        df = self.df_processed
        print(
            f" [{self._elapsed()}s] Starting preprocessing with pruning budget: {pruning_budget_k} kEUR..."
        )

        # 1. Prepare Data Arrays (for speed)
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

        # 2. Build Physical Graph (NetworkX)
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
                    # Only add edge once (u < v) to avoid duplicates
                    if u_idx < v_idx:
                        # Cost = Friction * Real Distance
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

        # 3. Calculate Pruning Cutoff
        # Logic: Max(15% of budget, or cost of ~6 median segments)
        median_cost = np.median(all_edge_costs) if all_edge_costs else 1
        budget_int = int(pruning_budget_k * self.SCALE_COST)
        limit_pct = budget_int * 0.15
        safety_floor = median_cost * 6
        cutoff_value = max(limit_pct, safety_floor)
        print(
            f"   [{self._elapsed()}s] Pruning Activated: Ignoring routes > {cutoff_value / self.SCALE_COST:.2f} kEUR"
        )

        # 4. Multi-Source Dijkstra (Path Calculation)
        print(f"   [{self._elapsed()}s] Calculating all viable routes (Dijkstra)...")
        paths_requirements = [[None for _ in range(n_species)] for _ in range(n_nodes)]

        for s_idx in range(n_species):
            # Identify source nodes (where species exists)
            sources = [i for i, val in enumerate(has_spec[s_idx]) if val]
            if not sources:
                continue

            # Compute shortest paths from ALL sources simultaneously
            # 'cutoff' prunes paths that are too expensive
            dists, paths = nx.multi_source_dijkstra(
                G, sources, weight="weight", cutoff=cutoff_value
            )

            for target, path_nodes in paths.items():
                if target in sources:
                    # Source cells don't need a path (cost 0)
                    paths_requirements[target][s_idx] = []
                    continue

                # Convert list of nodes [A, B, C] to list of edges [(A,B), (B,C)]
                edge_list = []
                for k in range(len(path_nodes) - 1):
                    u, v = sorted((path_nodes[k], path_nodes[k + 1]))
                    edge_list.append((u, v))
                paths_requirements[target][s_idx] = edge_list

        # Store pre-calculated data for the solver
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
        """
        Builds and solves the CP-SAT optimization model.

        Args:
            budget_k_euro (float): Total budget in thousands of euros.
            time_limit_sec (int): Max time for the solver to run (default 600s).

        Returns:
            pd.DataFrame: Solution dataframe with 'active_X', 'invest_X', and 'corridors' columns.
            None: If no feasible solution is found.
        """
        t_solve_start = time.time()

        # Ensure preprocessing is done
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

        # --- VARIABLES ---
        # x[i][s]: Cell i active for species s
        x = [
            [model.NewBoolVar(f"x_{i}_{s}") for s in range(n_species)]
            for i in range(n_nodes)
        ]

        # y[i][s]: Cell i receives investment for species s
        y = [
            [model.NewBoolVar(f"y_{i}_{s}") for s in range(n_species)]
            for i in range(n_nodes)
        ]

        # stress[i]: Indicator of predator-prey conflict in cell i
        stress = [model.NewBoolVar(f"strs_{i}") for i in range(n_nodes)]

        # Collect only relevant edges (those that are part of valid paths)
        used_edges = set()
        for i in range(n_nodes):
            for s in range(n_species):
                reqs = data["paths_requirements"][i][s]
                if reqs:
                    used_edges.update(reqs)

        # z[e]: Corridor edge e is built
        z = {}
        for e in used_edges:
            z[e] = model.NewBoolVar(f"z_{e}")

        # --- CONSTRAINTS ---
        obj_terms = []  # Objective function terms
        cost_terms = []  # Budget terms
        active_area_terms = [[] for _ in range(n_species)]  # For equity

        for i in range(n_nodes):
            # 1. Conflict: Martes and Eliomys cannot coexist
            model.Add(x[i][1] + x[i][2] <= 1)

            # 2. Stress: Martes + Oryctolagus => Stress
            model.Add(stress[i] >= x[i][1] + x[i][3] - 1)
            obj_terms.append(stress[i] * (-self.PENALTY_STRESS))

            for s in range(n_species):
                req_edges = data["paths_requirements"][i][s]

                # 3. Pruning: If path is None (too expensive), force x=0
                if req_edges is None:
                    model.Add(x[i][s] == 0)
                    continue

                # 4. Investment Logic
                model.Add(y[i][s] <= x[i][s])
                if not data["has_spec"][s][i]:
                    model.Add(x[i][s] <= y[i][s])  # Must invest if new
                cost_terms.append(y[i][s] * data["cost_adapt"][s][i])

                # 5. Path Implication: If x=1, all edges in path must be z=1
                if req_edges:
                    for e in req_edges:
                        model.AddImplication(x[i][s], z[e])

                # 6. Objective Terms
                # Benefit from active habitat
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

                # Benefit from restoration (Investment gain)
                added = self.RESTORED_Q - data["suitability"][s][i]
                if added > 0:
                    coef_y = int(
                        (added * self.W_VALS[s] * data["areas"][i] * self.SCALE_SCORE)
                        / self.SCALE_AREA
                    )
                    # Apply gain minus a small epsilon penalty (-1) to discourage trivial investments
                    obj_terms.append(y[i][s] * (coef_y - 1))
                else:
                    obj_terms.append(y[i][s] * -1)

                active_area_terms[s].append(x[i][s] * data["areas"][i])

        # 7. Corridor Costs
        for e, z_var in z.items():
            cost_terms.append(z_var * data["edges_unique"][e])
            obj_terms.append(z_var * -1)  # Small penalty for parsimony

        # 8. Budget Constraint (Safe against overflow for infinite budget)
        if budget_k_euro < 100_000:
            model.Add(sum(cost_terms) <= budget_int)

        # 9. Equity Constraints (Min/Max Share)
        total_active = sum(sum(active_area_terms[s]) for s in range(n_species))
        for s in range(n_species):
            s_sum = sum(active_area_terms[s])
            model.Add(s_sum * 100 >= self.EQUITY_MIN[s] * total_active)
            model.Add(s_sum * 100 <= self.EQUITY_MAX[s] * total_active)

        # --- EXECUTION ---
        model.Maximize(sum(obj_terms))

        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = time_limit_sec
        solver.parameters.num_search_workers = 12
        solver.parameters.log_search_progress = False

        print(f"[{self._elapsed()}s] Sending to solver (Limit: {time_limit_sec}s)...")
        print("   (Updating progress in real-time...)")

        callback = ProgressPrinter(time.time())
        status = solver.Solve(model, callback)
        print("\n")

        t_total_solve = time.time() - t_solve_start

        # --- RESULTS ---
        status_name = solver.StatusName(status)
        obj_val = solver.ObjectiveValue() / self.SCALE_SCORE

        print("\n" + "=" * 40)
        print(f" TOTAL TIME: {t_total_solve:.2f} seconds")
        if status == cp_model.OPTIMAL:
            print("OPTIMAL SOLUTION PROVEN")
        elif status == cp_model.FEASIBLE:
            print("FEASIBLE SOLUTION FOUND")
        else:
            print(f" FAILURE: {status_name}")
            return None

        print(f" Final Score: {obj_val:.2f}")
        print("=" * 40 + "\n")

        self.solution = self._format_solution(solver, x, y, z, data)
        return self.solution

    def _format_solution(self, solver, x, y, z, data):
        """Helper to convert solver variables into a readable DataFrame."""
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
        """Calculates and prints the real-world cost breakdown of the solution."""
        if self.solution is None:
            print(" Warning: No solution loaded.")
            return

        print("COST AUDIT")
        print("-" * 50)

        df_sol = self.solution
        df_base = self.df_processed
        df_merged = pd.merge(
            df_sol, df_base, on="grid_id", how="left", suffixes=("", "_base")
        )

        # 1. Adaptation Costs
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

        # 2. Corridor Costs (Geometric)
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
