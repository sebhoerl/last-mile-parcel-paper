import json

# This script writes out the configuration of an instance to solve.

if __name__ == "__main__":
    configuration = {}
    output_path = ""

    if "snakemake" in locals():
        output_path = snakemake.output[0]
        configuration = snakemake.params

        with open(output_path, "w+") as f:
            json.dump(configuration, f)