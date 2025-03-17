# Empirical data
#snakemake --unlock
#snakemake -j 10 -k --rerun-incomplete
# Simulated data
for FILE in config/simulated_mammals_te*.yaml;
do
  echo "Running ${FILE}"
  cp -f "${FILE}" config/config.yaml
  snakemake --unlock
  snakemake -s workflow/Simulations.smk -j 10 -k --rerun-incomplete
done
# Copy figures
./copy_artworks.sh