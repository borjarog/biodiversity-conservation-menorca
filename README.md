# Biodiversity Conservation Project in Menorca

## Project Description

This project implements an optimization model for biodiversity conservation on the island of Menorca. The objective is to determine the optimal allocation of resources (limited budget) to expand habitats for four key species through:

- **Cell adaptation**: Modifying cell habitats to make them suitable for specific species
- **Corridor construction**: Connecting existing populations with new areas through ecological corridors

### Target Species

1. *Atelerix algirus*
2. *Martes martes*
3. *Eliomys quercinus*
4. *Oryctolagus cuniculus*

## Project Structure
## Project Structure
```
biodiversity-conservation-menorca/
├── 1_documentation/                     # Model documentation and academic references
│   ├── model-design/                    # Model design versions (V1, V2, V3)
│   └── minimum-area-targets/            # Species-specific minimum area justifications
│       └── papers/                      # Scientific sources for area targets
│ 
├── 2_data/                              # Project datasets
│   ├── raw/                             # Original geographic data (dataset.geojson)
│   └── processed/                       # Processed datasets (CSV + GeoJSON)
│
├── 3_notebooks/                         # Jupyter notebooks (analysis & visualization)
│   ├── 00_project_statement.ipynb
│   ├── 01_eda_analysis.ipynb
│   ├── 02_visualization_maps.ipynb
│   ├── iterations/                      # Model iteration notebooks (V1/V2)
│   └── experiments/                     # Solver performance experiments
│       ├── datasets/                    # Benchmark results (multiple solvers)
│       ├── scripts/                     # Scripts used to run experiments
│       ├── analysis.ipynb               # Experiment analysis notebook
│       └── experiment visualizations    # PNG plots (comparisons, boxplots, etc.)
│
├── 4_source_code/                       # Core Python source code
│   ├── data_preparation.py              # Builds final dataset
│   ├── optimization_model.py            # CP-SAT optimization model
│   └── main.py                          # Model execution entry point
│
├── 5_results/                           # All generated outputs
│   ├── EDA/                             # Exploratory data analysis figures
│   ├── iterations/                      # Maps for V1 and V2 model runs
│   └── solutions/                       # Final optimal solution + species expansion maps
│       ├── species_expansion/
│       ├── optimal_solution_map.png
│       ├── solution_optimal.csv
│       ├── solution_v1.csv
│       └── solution_v2.csv
│
├── 6_final_submission/                  # Final formal delivery package
│
├── README.md
├── requirements.txt
├── pyproject.toml
└── uv.lock
```

## Model Evolution (Iterative Process)

The solution was developed through three distinct modeling phases to overcome scalability challenges:

### Phase 1: Initial Network Flow (V.1)
- **Concept:** A monolithic Mixed-Integer Linear Programming (MIP) model using flow variables ($f_{u,v,s}$) to enforce connectivity.
- **Outcome:** Computationally intractable on the full island grid due to the explosion of integer variables.

### Phase 2: Path Selection Strategy (V.2)
- **Concept:** Replaced flow variables with a pre-calculated path strategy. Connectivity was enforced via logical implications ($x \implies z$).
- **Outcome:** Solvable, but lacked geographical realism (corridors followed abstract grid adjacency rather than real terrain).

### Phase 3: Hybrid Geodetic Model (Final V.3)
- **Concept:** Integrates **real-world geography** (Haversine distances) and **graph pruning** in a preprocessing step using Dijkstra. The resulting optimized graph is fed into a compact **CP-SAT** model.
- **Outcome:** High performance, geographically accurate, and biologically robust solutions solved in seconds.

## Quick Start

### Prerequisites

Install all dependencies using the provided requirements file:

```bash
pip install -r requirements.txt
```

##  Model Execution & Workflow

The data has already been pre-processed and is ready to use. To run the optimization model:

```bash
cd 4_source_code
python main.py
```

Default execution settings:

- **Budget:** 1000 kEUR  
- **Time limit:** 600 seconds (10 minutes)

### Results

- The solution is saved at:  
  `5_results/solutions/solution_optimal.csv`
- A detailed **cost audit** is printed to the console.

---

## Workflow

### **1. Data Preparation (`data_preparation.py`)**

- Loads geographic data from:  
  `2_data/raw/dataset.geojson`
- Computes **suitability scores** based on land-cover type.
- Exports the final dataset in **CSV** and **GeoJSON** formats.

### **2. Optimization Model (`optimization_model.py`)**

- **Preprocessing:**  
  - Builds a georeferenced connectivity graph.  
  - Computes shortest paths (Dijkstra) from existing populations using real friction costs.

- **Optimization:**  
  - Uses **OR-Tools CP-SAT**.  
  - Ensures budget feasibility, biological equity, and species conflict management.

### **3. Analysis (`3_notebooks/`)**

- Notebooks for detailed exploration of model performance and resulting habitat networks.

---

## Model Configuration

### **Main Parameters (in `optimization_model.py`)**

- **Budget:** defined in `main.py` (default 1000 kEUR)  
- **Biological equity:** min and max active area per species  
- **Stress penalty:** penalizes species conflicts (Martes + Oryctolagus)

### **Model Constraints**

- **Budget:**  
  Total costs (adaptation + corridors) ≤ budget

- **Equity:**  
  Each species must occupy between X% and Y% of the total active area

- **Connectivity:**  
  Active cells must be connected to existing populations through paid corridors

- **Conflicts:**  
  Martes and Eliomys cannot coexist in the same cell

---

## Results

Outputs include:

- **Optimal solution:** CSV with selected cells and active species  
- **Cost audit:** expense breakdown (adaptation vs. corridors)  
- **Performance metrics:** objective value, optimality gap, runtime

## Optimal Solution Map

![Optimal Solution Map: Active Habitats and Ecological Corridors](5_results/solutions/optimal_solution_map.png)

---

## Additional Documentation

For more technical details:

- `1_documentation/model-design/`: full PDF documentation of model evolution  
- `3_notebooks/01_eda_analysis.ipynb`: exploratory analysis notebook

---

## Author

Academic project on biodiversity conservation in Menorca.

