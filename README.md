
**Compiled binaries and instructions for BayesCode are available
at [github.com/ThibaultLatrille/bayescode](https://github.com/ThibaultLatrille/bayescode)**

# ChronoPhylogram

This repository is meant to provide the necessary scripts and data to reproduce the figures shown in the manuscript.
The experiments can either run on a local computer or in a cluster configuration (slurm).

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

Install [BayesCode](https://github.com/ThibaultLatrille/bayescode) to obtain the executable `nodetraits`
and `readnodetraits`:

```
conda install -c bioconda -c conda-forge bayescode
```

Alternatively, you can compile BayesCode from source in the `utils` folder:

```bash
mkdir utils
cd utils
git clone https://github.com/ThibaultLatrille/BayesCode
cd BayesCode
make tiny
``` 

## 2. Run empirical analysis

Run `snakemake`:

```
snakemake -j 8 -k
```
