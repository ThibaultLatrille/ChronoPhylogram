for FILE in config/simulated_*.yaml;
do
  echo "Running ${FILE}"
  cp -f "${FILE}" config/config.yaml
  snakemake --unlock
  snakemake -s workflow/Simulations.smk -j 10 -k --rerun-incomplete
done
