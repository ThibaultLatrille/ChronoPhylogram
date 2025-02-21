import os
import argparse
import numpy as np
import pandas as pd
from os.path import isdir
from libraries_plot import vert_boxplot


def format_label(l):
    l = l.replace("b ", "").replace("m ", "♂").replace("f ", "♀ ")
    l = l.replace("bodyMass", "body mass").replace("brainMass", "brain mass")
    return l.replace("REML", "(REML)").replace("nodes", "")


def main(folder, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    models = {"REML": "simple_BM_SwitchREML", "nodes": "simple_BM_Switchnodes"}
    x_input = {}
    for directory in sorted(os.listdir(folder)):
        if not isdir(f"{folder}/{directory}"):
            continue
        for rb_name, rb in models.items():
            path = f"{folder}/{directory}/{rb}.log.gz"
            if not os.path.exists(path):
                print(path, "not found")
                continue
            label = f"{directory}_{rb_name}"
            p_nuc = pd.read_csv(path, sep="\t")["is_nuc"]
            # cut in 10 slices
            p_slice = np.array_split(p_nuc, 10)
            print(len(p_nuc), "points in", label)
            # compute the mean across 100 points
            x_input[label] = [np.mean(p) for p in p_slice]

    vert_boxplot(x_input, "Probability of Phylogram", output, yscale="uniform",
                 format_label=format_label, empirical=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
