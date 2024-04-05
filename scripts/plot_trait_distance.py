import argparse
from collections import defaultdict
import numpy as np
from libraries import *
import matplotlib.pyplot as plt


def main(input_distance_tree: str, input_simu_tree_list: list[str], output_path: str):
    d_tree = open_tree(input_distance_tree, format_ete3=1)
    dict_x, dict_y = defaultdict(list), defaultdict(list)
    for input_simu_tree in input_simu_tree_list:
        simu_tree = open_tree(input_simu_tree, format_ete3=1)
        for n_t, n_d in zip(simu_tree.traverse(), d_tree.traverse()):
            if n_t.is_root() or n_t.up.is_root():
                continue
            assert n_t.name == n_d.name, f"Taxon names do not match: {n_t.name} != {n_d.name}"
            dict_x[n_t.name].append(np.sqrt(n_d.dist))
            trait = float(getattr(n_t, "Phenotype_mean"))
            trait_parent = float(getattr(n_t.up, "Phenotype_mean"))
            dict_y[n_t.name].append(np.abs(trait - trait_parent))
    x = np.array([np.mean(v) for v in dict_x.values()])
    x_std = np.array([np.std(v) for v in dict_x.values()])
    y = np.array([np.mean(v) for v in dict_y.values()])
    y_std = np.array([np.std(v) for v in dict_y.values()])
    assert len(x) == len(y), "Different number of x and y values"
    # Scatter plot
    plt.figure(figsize=(12, 8))
    plt.scatter(x, y, alpha=0.4)
    # Display the range of the inference nodes
    plt.errorbar(x, y, xerr=x_std, yerr=y_std, fmt='o', color='black', alpha=0.5)
    # regression line and r-squared
    a, b = np.polyfit(x, y, 1)
    r2 = np.corrcoef(x, y)[0, 1] ** 2
    x_lim = np.array([min(0, min(x)), max(x)])
    plt.plot(x_lim, a * x_lim + b, color="red", label=f"y = {a:.2f}x + {b:.2f} (R² = {r2:.2f})")
    # Horizontal line at 0
    plt.axhline(0, color="black", linestyle="--")
    plt.xlabel("sqrt(d)")
    plt.ylabel("|Dz|")
    plt.legend()
    plt.savefig(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--distance_tree", help="Input distance tree file", required=True)
    parser.add_argument("--simu_tree", help="Input simulation tree file", nargs="+", required=True)
    parser.add_argument("--output", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.distance_tree, args.simu_tree, args.output)
