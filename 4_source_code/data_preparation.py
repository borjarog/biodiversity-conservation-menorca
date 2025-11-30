import os
import sys
import geopandas as gpd

# ----------------------------------------------------
# PATHS
# ----------------------------------------------------

input_data_path = os.path.join(
    os.path.dirname(__file__), "../2_data/raw/dataset.geojson"
)
output_data_path = os.path.join(
    os.path.dirname(__file__), "../2_data/processed/final_dataset.geojson"
)

# ----------------------------------------------------
# SUITABILITY SCORE MAP
# ----------------------------------------------------

score_map = {
    "High": 3.0,
    "Moderate-High": 2.5,
    "Moderate": 2.0,
    "Low-Moderate": 1.5,
    "Low": 1.0,
    "Very Low": 0.0,
}

# Species rules
from_rules = {
    "atelerix": {
        "Discontinuous Urban Fabric": "Moderate",
        "Continuous Urban Fabric": "Low",
        "Industrial or Commercial Units": "Very Low",
        "Airports": "Very Low",
        "Port Areas": "Very Low",
        "Sport and Leisure Facilities": "Very Low",
        "Pastures": "High",
        "Non-irrigated Arable Land": "High",
        "Permanently Irrigated Land": "Low",
        "Complex Cultivation Patterns": "High",
        "Land Principally Occupied by Agriculture with Significant Areas of Natural Vegetation": "High",
        "Sclerophyllous Vegetation": "High",
        "Transitional Woodland-Shrub": "High",
        "Natural Grasslands": "Moderate-High",
        "Broad-leaved Forests": "Moderate",
        "Mixed Forests": "Moderate",
        "Coniferous Forests": "Low",
        "Peatbogs": "Very Low",
        "Inland Marshes": "Very Low",
        "Coastal Lagoons": "Very Low",
        "Estuaries": "Very Low",
        "Intertidal Flats": "Very Low",
        "Water Courses": "Low-Moderate",
    },
    "martes": {
        "Discontinuous Urban Fabric": "Low",
        "Continuous Urban Fabric": "Very Low",
        "Industrial or Commercial Units": "Very Low",
        "Airports": "Very Low",
        "Port Areas": "Very Low",
        "Sport and Leisure Facilities": "Very Low",
        "Pastures": "Low",
        "Non-irrigated Arable Land": "Low",
        "Permanently Irrigated Land": "Very Low",
        "Complex Cultivation Patterns": "Moderate",
        "Land Principally Occupied by Agriculture with Significant Areas of Natural Vegetation": "Moderate-High",
        "Sclerophyllous Vegetation": "High",
        "Transitional Woodland-Shrub": "High",
        "Natural Grasslands": "Low-Moderate",
        "Broad-leaved Forests": "High",
        "Mixed Forests": "High",
        "Coniferous Forests": "Moderate-High",
        "Peatbogs": "Very Low",
        "Inland Marshes": "Very Low",
        "Coastal Lagoons": "Very Low",
        "Estuaries": "Very Low",
        "Intertidal Flats": "Very Low",
        "Water Courses": "Low-Moderate",
    },
    "eliomys": {
        "Discontinuous Urban Fabric": "Low",
        "Continuous Urban Fabric": "Very Low",
        "Industrial or Commercial Units": "Very Low",
        "Airports": "Very Low",
        "Port Areas": "Very Low",
        "Sport and Leisure Facilities": "Very Low",
        "Pastures": "Low",
        "Non-irrigated Arable Land": "Low",
        "Permanently Irrigated Land": "Very Low",
        "Complex Cultivation Patterns": "Moderate",
        "Land Principally Occupied by Agriculture with Significant Areas of Natural Vegetation": "Moderate-High",
        "Sclerophyllous Vegetation": "High",
        "Transitional Woodland-Shrub": "High",
        "Natural Grasslands": "Low",
        "Broad-leaved Forests": "High",
        "Mixed Forests": "High",
        "Coniferous Forests": "Moderate-High",
        "Peatbogs": "Very Low",
        "Inland Marshes": "Very Low",
        "Coastal Lagoons": "Very Low",
        "Estuaries": "Very Low",
        "Intertidal Flats": "Very Low",
        "Water Courses": "Low-Moderate",
    },
    "oryctolagus": {
        "Discontinuous Urban Fabric": "Low",
        "Continuous Urban Fabric": "Very Low",
        "Industrial or Commercial Units": "Very Low",
        "Airports": "Very Low",
        "Port Areas": "Very Low",
        "Sport and Leisure Facilities": "Very Low",
        "Pastures": "High",
        "Non-irrigated Arable Land": "High",
        "Permanently Irrigated Land": "Low",
        "Complex Cultivation Patterns": "High",
        "Land Principally Occupied by Agriculture with Significant Areas of Natural Vegetation": "High",
        "Sclerophyllous Vegetation": "High",
        "Transitional Woodland-Shrub": "High",
        "Natural Grasslands": "High",
        "Broad-leaved Forests": "Moderate",
        "Mixed Forests": "Moderate",
        "Coniferous Forests": "Low-Moderate",
        "Peatbogs": "Very Low",
        "Inland Marshes": "Very Low",
        "Coastal Lagoons": "Very Low",
        "Estuaries": "Very Low",
        "Intertidal Flats": "Very Low",
        "Water Courses": "Low",
    },
}

# ----------------------------------------------------
# GET SCORE FUNCTION AND XY EXTRACTION
# ----------------------------------------------------


def get_score(land, species):
    label = from_rules[species].get(land, "Very Low")
    return score_map.get(label, 0.0)


def extract_xy(grid_id):
    try:
        _, x, y = grid_id.split("_")
        return int(x), int(y)
    except Exception:
        return None, None


# ----------------------------------------------------
# GET NEIGHBORS FUNCTION
# ----------------------------------------------------


def compute_neighbors(df):
    coord_to_id = {(row.grid_x, row.grid_y): row.grid_id for _, row in df.iterrows()}

    neighbor_list = []

    for _, row in df.iterrows():
        x, y = row.grid_x, row.grid_y

        candidates = [
            (x - 1, y),
            (x + 1, y),
            (x, y - 1),
            (x, y + 1),
            (x - 1, y - 1),
            (x - 1, y + 1),
            (x + 1, y - 1),
            (x + 1, y + 1),
        ]

        neigh_ids = []
        for c in candidates:
            if c in coord_to_id:
                neigh_ids.append(coord_to_id[c])

        neighbor_list.append(";".join(neigh_ids))

    df["neighbors"] = neighbor_list
    return df


# ----------------------------------------------------
# MAIN PIPELINE
# ----------------------------------------------------

if __name__ == "__main__":
    if not os.path.exists(input_data_path):
        print(f"[ERROR] File not found: {input_data_path}")
        sys.exit(1)

    gdf = gpd.read_file(input_data_path)
    print(f"[INFO] Loaded dataset: {len(gdf)} cells.")

    # Extract grid_x and grid_y
    gdf["grid_x"], gdf["grid_y"] = zip(*gdf["grid_id"].apply(extract_xy))
    gdf["grid_x"], gdf["grid_y"] = zip(*gdf["grid_id"].apply(extract_xy))

    gdf = compute_neighbors(gdf)

    print("[INFO] Extracted grid_x and grid_y from grid_id.")
    print("[INFO] Computed neighbors for each cell.")

    # Compute suitability columns
    for sp in ["atelerix", "martes", "eliomys", "oryctolagus"]:
        gdf[f"suitability_{sp}"] = gdf["dominant_land_cover_name"].apply(
            lambda x: get_score(x, sp)
        )

    # Save final dataset
    print(f"[INFO] Saving final GeoJSON: {output_data_path}")
    gdf.to_file(output_data_path, driver="GeoJSON")

    # Export also to CSV (geometry becomes WKT by default unless you drop it)
    csv_output_path = os.path.join(
        os.path.dirname(__file__), "../2_data/processed/final_dataset.csv"
    )

    print(f"[INFO] Saving CSV: {csv_output_path}")
    gdf.to_csv(csv_output_path, index=False)

    print("[INFO] DONE.")
