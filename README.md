# Biodiversity Conservation Project in Menorca

## 📋 Project Description

This project implements an optimization model for biodiversity conservation on the island of Menorca. The objective is to determine the optimal allocation of resources (limited budget) to expand habitats for four key species through:

- **Cell adaptation**: Modifying cell habitats to make them suitable for specific species
- **Corridor construction**: Connecting existing populations with new areas through ecological corridors

### Target Species

1. **Atelerix algirus** (Algerian hedgehog) - Weight: 1.0
2. **Martes martes** (Pine marten) - Weight: 1.0
3. **Eliomys quercinus** (Garden dormouse) - Weight: 2.0
4. **Oryctolagus cuniculus** (European rabbit) - Weight: 1.5

## 🏗️ Project Structure

```text
biodiversity-conservation-menorca/
├── 1_documentation/          # Model documentation and justifications
│   ├── model-design/         # Model design versions
│   └── minimum-area-targets/ # Minimum area target justifications
├── 2_data/                   # Project data
│   ├── raw/                  # Original data (dataset.geojson)
│   └── processed/            # Processed data (final_dataset.csv/geojson)
├── 3_notebooks/              # Exploratory analysis and visualizations
│   ├── 00_project_statement.ipynb.ipynb
│   ├── 01_eda_analysis.ipynb
│   ├── 02_visualization_maps.ipynb
│   └── iterations/           # Model iterations (v1, v2)
├── 4_source_code/            # Main source code
│   ├── data_preparation.py   # Data preparation and processing
│   ├── optimization_model.py # CP-SAT optimization model
│   └── main.py               # Main execution script
├── 5_results/                # Generated results
│   ├── EDA/                  # Exploratory data analysis
│   ├── solutions/            # Optimal solutions and maps
│   ├── iterations/           # Iteration results
│   └── tables/               # Result tables
└── 6_final_submission/       # Final submission
```

## 🚀 Quick Start

### Prerequisites

```bash
pip install pandas geopandas numpy networkx ortools matplotlib seaborn
```

Or install from `requirements.txt` (if available).

### Model Execution

1. **Prepare the data** (if not already processed):

```bash
cd 4_source_code
python data_preparation.py
```

2. **Run the optimization model**:

```bash
python main.py
```

   **Note**: The script will run the solver with a budget of 1000 kEUR and a time limit of 600 seconds (10 minutes).

3. **Results**:

   The solution is saved to `5_results/solutions/solution_optimal.csv` and a cost audit is generated in the console.

## 📊 Workflow

1. **Data Preparation** (`data_preparation.py`):
   - Loads geographic data from `2_data/raw/dataset.geojson`
   - Calculates suitability scores for each species based on land cover type
   - Identifies neighbors for each cell
   - Exports final dataset in CSV and GeoJSON formats

2. **Optimization Model** (`optimization_model.py`):
   - Builds connectivity graph between cells
   - Calculates shortest paths (Dijkstra) from existing populations
   - Solves optimization problem with CP-SAT (Google OR-Tools)
   - Considers budget constraints, biological equity, and species conflicts

3. **Analysis and Visualization** (`3_notebooks/`):
   - Exploratory data analysis (EDA)
   - Map and result visualizations
   - Model iteration comparisons

## 🔧 Model Configuration

### Main Parameters (in `optimization_model.py`)

- **Budget**: Defined in `main.py` (default: 1000 kEUR)
- **Species weights**: `W_VALS = [1.0, 1.0, 2.0, 1.5]`
- **Biological equity**: Minimum and maximum ranges of active area per species
- **Stress penalty**: Penalizes species conflicts (Martes + Oryctolagus)

### Model Constraints

1. **Budget**: Total costs (adaptation + corridors) ≤ budget
2. **Equity**: Each species must have between X% and Y% of total active area
3. **Connectivity**: Active cells must be connected to existing populations
4. **Conflicts**: Martes and Eliomys cannot coexist in the same cell

## 📈 Results

Results include:

- **Optimal solution**: CSV with selected cells and active species
- **Visualization maps**: Habitat expansion by species
- **Cost audit**: Breakdown of expenses in adaptation and corridors
- **Performance metrics**: Objective score, optimality gap, execution time

## 📚 Additional Documentation

For more technical details, see:

- `DOCUMENTATION.md`: Complete technical documentation
- `1_documentation/model-design/`: Model design and evolution
- `3_notebooks/01_eda_analysis.ipynb`: Detailed exploratory analysis

## 👥 Author

Academic project on biodiversity conservation in Menorca.
