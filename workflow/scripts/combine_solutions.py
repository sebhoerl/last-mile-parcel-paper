import pandas as pd

# This script takes the individual outcomes (stops, legs, tours, ...) of multiple
# scenarios and combines them into aggregated data frames

if __name__ == "__main__":
    stops_paths = []
    legs_paths = []
    tours_paths = []
    vehicles_paths = []
    instances_paths = []

    for name in snakemake.input.keys():
        if name.startswith("stops:"):
            stops_paths.append(snakemake.input[name])

        if name.startswith("legs:"):
            legs_paths.append(snakemake.input[name])

        if name.startswith("tours:"):
            tours_paths.append(snakemake.input[name])

        if name.startswith("vehicles:"):
            vehicles_paths.append(snakemake.input[name])

        if name.startswith("instances:"):
            instances_paths.append(snakemake.input[name])
    
    df_stops = pd.concat([
        pd.read_parquet(path) for path in stops_paths
    ])

    df_legs = pd.concat([
        pd.read_parquet(path) for path in legs_paths
    ])

    df_tours = pd.concat([
        pd.read_parquet(path) for path in tours_paths
    ])

    df_vehicles = pd.concat([
        pd.read_parquet(path) for path in vehicles_paths
    ])

    df_instances = pd.concat([
        pd.read_parquet(path) for path in instances_paths
    ])

    df_stops.to_parquet(snakemake.output["stops"])
    df_legs.to_parquet(snakemake.output["legs"])
    df_tours.to_parquet(snakemake.output["tours"])
    df_vehicles.to_parquet(snakemake.output["vehicles"])
    df_instances.to_parquet(snakemake.output["instances"])
