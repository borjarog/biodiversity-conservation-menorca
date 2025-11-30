"""
Menorca Conservation Project - Main Execution Script
====================================================

This script serves as the primary entry point for the biodiversity conservation
optimization pipeline. It orchestrates the entire process:
1. Locates and verifies input data (CSV and GeoJSON).
2. Initializes the optimization solver (MenorcaSolver).
3. Executes the optimization with defined budget and time limits.
4. Audits the results and saves the optimal solution to disk.

Usage:
    Run this script directly from the terminal:
    $ python main.py

Directory Structure Assumption:
    This script assumes a specific project structure where data resides in '../2_data'
    and results are saved to '../5_results'.
"""

import os
import sys
import traceback

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from optimization_model import MenorcaSolver


def main():
    """
    Orchestrates the optimization workflow.

    Steps:
        1. Sets up file paths for input data and output results.
        2. Validates the existence of required data files.
        3. Instantiates the MenorcaSolver.
        4. Runs the solver with a budget of 1000 k€ and a 10-minute time limit.
        5. Prints a cost audit and saves the solution to CSV if successful.
    """
    print("STARTING SOLVER EXECUTION")
    print("================================")

    # 1. FILE PATH SETUP
    # -------------------------------------------------------
    # Calculate the base directory (project root) relative to this script
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    csv_path = os.path.join(base_dir, "2_data", "processed", "final_dataset.csv")
    geo_path = os.path.join(base_dir, "2_data", "processed", "final_dataset.geojson")

    # Verify inputs before starting expensive operations
    if not os.path.exists(csv_path) or not os.path.exists(geo_path):
        print("ERROR: Data files not found.")
        print(f"Checked CSV path: {csv_path}")
        print(f"Checked GeoJSON path: {geo_path}")
        return

    try:
        # 2. INITIALIZATION
        # -------------------------------------------------------
        print(f"Initializing solver with data from: {base_dir}")
        solver = MenorcaSolver(csv_path, geo_path)

        # 3. EXECUTION CONFIGURATION
        # -------------------------------------------------------
        # Budget in thousands of euros (k€)
        BUDGET_K_EUR = 1000
        # Max execution time in seconds (600s = 10 minutes)
        TIME_LIMIT_SEC = 600

        # Run the optimization
        df_result = solver.solve(
            budget_k_euro=BUDGET_K_EUR, time_limit_sec=TIME_LIMIT_SEC
        )

        # 4. RESULTS HANDLING
        # -------------------------------------------------------
        if df_result is not None:
            # Print detailed cost breakdown to console
            solver.audit_results()

            # Define output directory
            results_dir = os.path.join(base_dir, "5_results", "solutions")

            # Create directory if it doesn't exist (Safety check)
            if not os.path.exists(results_dir):
                os.makedirs(results_dir)
                print(f"\nCreated output directory: {results_dir}")

            # Save to CSV
            out_file = os.path.join(results_dir, "solution_optimal.csv")
            df_result.to_csv(out_file, index=False)
            print(f"\nSolution saved successfully: {out_file}")
        else:
            print("\nNo optimal or feasible solution was found within the time limit.")

    except Exception as e:
        print("\nERROR DURING EXECUTION:")
        print(str(e))
        traceback.print_exc()


if __name__ == "__main__":
    main()
