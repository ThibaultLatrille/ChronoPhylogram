#!/usr/bin/env python3
import os
import argparse
from gzip import open as gzopen
import pandas as pd
from ete3 import Tree
from Bio import Phylo
from libraries import open_tree, prune_tree, replace_last


def main(input_tree: str, input_traits: str, output_tree: str, output_traits: str):
    for path in [input_traits, input_tree]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [output_traits, output_tree]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    # Convert newick to nexus including root
    Phylo.convert(input_tree, "newick", output_tree, "nexus")
    traits = pd.read_csv(input_traits, sep="\t")
    with open(output_traits, "w") as f:
        f.write('#NEXUS\n\n')
        f.write("Begin data;\n")
        f.write(f"Dimensions ntax={len(traits)} nchar=1;\n")
        f.write("Format datatype=Continuous missing=? gap=-;\n")
        f.write("Matrix\n")
        for row in traits.itertuples():
            f.write(f"{row.TaxonName}\t{row.Phenotype_mean}\n")
        f.write(";\n")
        f.write("End;\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_tree", help="Input tree file", required=True)
    parser.add_argument("--input_traits", help="Input traits file", required=True)
    parser.add_argument("--output_tree", help="Output tree file", required=True)
    parser.add_argument("--output_traits", help="Output traits file", required=True)
    args = parser.parse_args()
    main(args.input_tree, args.input_traits, args.output_tree, args.output_traits)
