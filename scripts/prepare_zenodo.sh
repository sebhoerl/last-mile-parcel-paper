set -e

if [ ! -f README.md ]; then
  echo "This script should be executed in the root directory for the repository"
  exit 1
fi

# create directory
rm -rf zenodo
mkdir -p zenodo

cp results/depots/prepared.gpkg zenodo

for year in 2019 2024 2030; do
  cp resources/external/population_${year}/lyon_${year}_100pct_homes.gpkg zenodo
  cp resources/external/population_${year}/lyon_${year}_100pct_persons.csv zenodo
done
