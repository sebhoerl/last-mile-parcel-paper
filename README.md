# Decarbonization policies for last-mile parcel delivery

This repository contains a open-data-driven and replicable for the modeling of last mile parcel deliveries in a city. The model is presented for the use case of Lyon with fully reproducible results.

The model structure and assumptions are described in the following research article:

> Hörl, S. Sebastian Hörl, Jakob Puchinger. Decarbonization policies for last-mile parcel delivery: An adaptable open-data case study for Lyon.

While most information on the individual modeling steps can be found in the article, here instructions on how to run the model and reproduce the results shall be given.

The model is structured as a consistent pipeline from raw data to final outputs based on [Snakemake](https://snakemake.readthedocs.io/en/stable/). As soon as the input data is put in the right place, the model should be able to run without further intervention. However, solving one scenario takes about 2-3 days on a standard machine, so running all scenarios requires about two weeks. Further information on the role of each notebook and script is given in the respective files.

## Running the model

To run the model, a couple of steps need to be followed. The following sections explain how to set up the runtime environment, install the Python dependencies, collect the input data sets, execute the pipeline, and run the analysis code.

### Linux environment

It is recommended to run the model on a Linux machine or in the Linux subsystem for Windows (WSL). The main reason is that the Vehicle Routing Problem (VRP) solver [VROOM](https://github.com/VROOM-Project/vroom) is used heavily inside the pipeline. However, we require a customized version which will be pulled automatically from Github and built using `cmake`.

To set up the common dependencies that are required for building VROOM on a standard Ubuntu image, have a look at `environment.sh`.

Additionally, `osmosis` needs to be accessible from the command line to filter large OpenStreetMap data sets. A recent version of the tool can be downloaded [here](https://github.com/openstreetmap/osmosis/releases/tag/0.49.2).

### Python environment

The whole pipeline is based on Python with various package dependencies. This repository contains `environment.yml` which describes a `conda` / `mamba` environment. To set up the environment, run

```bash
mamba env create -f environment.yml
```

A new environment called `parcels` is set up inside which the model pipeline must be called. This means that whenever you run the following pipeline code, the environment needs to be activated using

```bash
mamba activate parcels
```

### Data collection


To run the model, first some raw data needs to be collected and placed into the `resources` directoy. The following sections describe how to obtain this data.

#### OpenSteetMap

To obtain a road network of the region, a cut-out from OpenStreetMap is required. Such a cut-out can be obtained from [Geofabrik](https://download.geofabrik.de/europe/france/rhone-alpes.html) for the Rhône-Alpes region around Lyon. To make sure the data is the same as used in the accompanying paper, go to the [cut-out archive](https://download.geofabrik.de/europe/france/) and download the file with the name `rhone-alpes-220101.osm.pbf` from January 2022. Place the file into `resources/external/osm`.

#### IRIS zoning system

Furthermore, the geographic description of the zoning system used by the French statistical office (INSEE) is required. It can be downloaded from [IGN](https://geoservices.ign.fr/contoursiris). After navigating to the website, download the 2022 edition and place the file `CONTOURS-IRIS_2-1__SHP__FRA_2022-01-01.7z` into `resources/external/iris`.

#### SIRENE enterprise information

In case you want to run the pipeline from start to end, you will need to download the French enterprise registry SIRENE. It can be downloaded from [INSEE](https://www.data.gouv.fr/fr/datasets/base-sirene-des-entreprises-et-de-leurs-etablissements-siren-siret/). Download *Sirene: Ficher StockEtablissement du 01 [Month Year]* and put the file `StockEtablissement_utf8.zip` into `resources/external/sirene`.

The rule `prepare_depots` inside the workflow will process this file and use a public API accessing the French national address database (BAN) to geolocate the contained addresses. However, it should be noted that the SIRENE data linked above is updated regularly so it is not possible to use the exact same input as was used in the research paper.

**Alternatively**, to maintain reproducibility, you can download the output of the `prepare_depots` rule from [Zenodo](https://doi.org/10.5281/zenodo.13933045). Download the file `prepared.gpkg` and place it at `results/depots/prepared.gpkg`. This way, execution of the `prepare_depots` rule will not be triggered when executing the pipeline, but the precalculated information from this file will be used.

#### Synthetic population

Finally, synthetic population data is needed. The synthetic population for Lyon is generated using the replicable method described in 

> Hörl, S. and M. Balac (2021) [Synthetic population and travel demand for Paris and Île-de-France based on open and publicly available data](https://www.sciencedirect.com/science/article/pii/S0968090X21003016), Transportation Research Part C, 130, 103291.

Applying that process requires itself the collection of about twelve individual open data sets in France and then running a pipeline for multiple hours. The whole process is [documented in detail](https://github.com/eqasim-org/ile-de-france/blob/develop/docs/population.md) and can be followed to generate the input data for the present parcel delivery workflow. Of interest are the files `homes.gpkg` and `persons.csv`. The [process for Lyon](https://github.com/eqasim-org/ile-de-france/blob/develop/docs/cases/lyon.md) should be followed and population projections [should be activated](https://github.com/eqasim-org/ile-de-france/blob/develop/docs/population.md#population-projections). The resulting three populations for Lyon should then be put into `resources` using the following names by replacing `[year]` with the three projection years 2019, 2024, and 2030:

- `external/population_[year]/lyon_[year]_100pct_homes.gpkg`
- `external/population_[year]/lyon_[year]_100pct_persons.csv`

**Alternatively**, you may download pregenerated synthetic populations from [Zenodo](https://doi.org/10.5281/zenodo.13933045). Copy the prefixed `homes.gpkg` and `persons.csv` files into `external` to arrive at the directory structure outlined above.

### Pipeline execution

To run the model, activate the Python environment and call the modeling pipeline using `snakemake` in the root folder of this repository. It will pick up the *rules* to generate all the intermediate and final files from `workflow/Snakefile`:

```bash
snakemake -c 12 all
```

You can define the number of cores to use by providing a different number after `-c`. If all input data has been put in the right place, the pipeline should start without any error. Note that to run the pipeline an Internet connection is necessary. In particular, the pipeline will clone `vroom` from a Github repository and build it (`build_vroom` rule), and it will use the French address API to geocode the locations of distribution centers (`prepare_depots` rule).

The final results for each scenario will be placed into `results/scenario_solutions`. The folder contains analysis files on different levels of aggregation: per instance (depot), per vehicle leg, per stop, per tour, per vehicle.

### Analysis

The output can be analyzed by running the notebooks in `paper`. They will display the results directly and generate the figures and tables that are used in the corresponding paper in `paper/figures` and `paper/tables`, respectively.

To generate the maps, the `Prepare Map` notebook generates several geographic files that can be visualized using [QGIS](https://qgis.org/). A QGIS project file is given in `paper/map/paper.qgz`.

## Acknowledgements

This model has been partly developed in the scope of the project LEAD, which has received funding from the European Union's Horizon 2020 research and innovation program under grant agreement no. 861598. The content of this documentation does not reflect the official opinion of the European Union. Responsibility for the information and views expressed in this repository lies entirely with the authors.
