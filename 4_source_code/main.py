import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from optimization_model import MenorcaSolver


def main():
    print("STARTING SOLVER EXECUTION")
    print("================================")

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    csv_path = os.path.join(base_dir, "2_data", "processed", "final_dataset.csv")
    geo_path = os.path.join(base_dir, "2_data", "processed", "final_dataset.geojson")

    if not os.path.exists(csv_path) or not os.path.exists(geo_path):
        print("ERROR: Data files not found.")
        print(f"Checked path: {csv_path}")
        return

    try:
        # We instantiate the Solver
        solver = MenorcaSolver(csv_path, geo_path)

        # We define budget and time limit
        budget = 1000
        time_limit = 600

        df_result = solver.solve(budget_k_euro=budget, time_limit_sec=time_limit)

        # If a solution was found, we proceed to save it
        if df_result is not None:
            solver.audit_results()
            results_dir = os.path.join(base_dir, "5_results", "solutions")

            if not os.path.exists(results_dir):
                os.makedirs(results_dir)
                print(f"\nCreated directory: {results_dir}")

            out_file = os.path.join(results_dir, "solution_optimal.csv")
            df_result.to_csv(out_file, index=False)
            print(f"\nSolution saved successfully: {out_file}")

    except Exception as e:
        print(f"\nERROR DURING EXECUTION:")
        print(e)
        import traceback

        traceback.print_exc()


if __name__ == "__main__":
    main()
