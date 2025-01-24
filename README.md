# ChronoPhylogram

This repository is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript.
The experiments are meant to run on Linux/Unix/MacOS operating systems.

If problems and/or questions are encountered, feel free
to [open issues](https://github.com/ThibaultLatrille/ChronoPhylogram/issues).

## 0. Local copy

Clone the repository and `cd` to the dir.

```
git clone https://github.com/ThibaultLatrille/ChronoPhylogram
cd ChronoPhylogram
```

## 1. Installation

### General dependencies

Install python3 packages

```
sudo apt install -qq -y python3-dev python3-pip
pip3 install snakemake scipy numpy matplotlib pandas ete3 bio statsmodels --user
```

Install [RevBayes](https://revbayes.github.io/) to obtain the executable `rb`:
With conda:
```
conda install -c bioconda -c conda-forge revbayes
```

## 2. Run empirical analysis

Run `snakemake`:

```
snakemake -j 8 -k
```

## 2. Run simulated analysis

The simulation program is necessary to generate the data for the simulated analysis, which are stored in the `data_simulated` folder.\
The details to compile the simulation program and its usage are described in the `README.md` file in the `utils/simulator` folder.

Run `snakemake`:

```
snakemake -s workflow/Simulations.smk -j 8 -k 
```

