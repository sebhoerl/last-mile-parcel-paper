set -e

git clone ...
cd ..
mvn clean package -Pstandalone --also-make -DskipTests=true --projects ile_de_france
