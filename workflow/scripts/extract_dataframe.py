import geopandas as gpd

# This script extracts a smaller dataframe from a larger one by filtering on one column

def extract_dataframe(source_path, target_path, column):
    df_source = gpd.read_file(source_path)
    df_target = df_source[df_source[column].astype(str) == str(value)]
    df_target.to_file(target_path)

if __name__ == "__main__":
    source_path = "results/parcels/baseline_2022/assigned.gpkg"
    target_path = "results/parcels/baseline_2022/operators/chronopost.gpkg"
    column = "operator"
    value = "chronopost"

    if "snakemake" in locals():
        source_path = snakemake.input[0]
        target_path = snakemake.output[0]
        column = snakemake.params["column"]
        value = snakemake.params["value"]

    extract_dataframe(source_path, target_path, column)
