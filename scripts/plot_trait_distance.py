import argparse
import numpy as np
from libraries import *
import matplotlib.pyplot as plt


def main(input_distance_tree: str, input_simu_tree: str, output_path: str):
    d_tree = open_tree(input_distance_tree, format_ete3=1)
    simu_tree = open_tree(input_simu_tree, format_ete3=1)
    x, y = [], []
    for n_t, n_d in zip(simu_tree.traverse(), d_tree.traverse()):
        if n_t.is_root() or n_t.up.is_root():
            continue
        assert n_t.name == n_d.name, f"Taxon names do not match: {n_t.name} != {n_d.name}"
        x.append(np.sqrt(n_d.dist))
        trait = float(getattr(n_t, "Phenotype_mean"))
        trait_parent = float(getattr(n_t.up, "Phenotype_mean"))
        y.append(np.abs(trait - trait_parent))
    assert len(x) == len(y), "Different number of x and y values"
    # Scatter plot
    plt.figure(figsize=(12, 8))
    plt.scatter(x, y)
    # regression line and r-squared
    m, b = np.polyfit(x, y, 1)
    r2 = np.corrcoef(x, y)[0, 1] ** 2
    plt.plot(x, m * np.array(x) + b, color="red", label=f"y = {m:.2f}x + {b:.2f} (R² = {r2:.2f})")
    plt.xlabel("sqrt(d)")
    plt.ylabel("|Dz|")
    plt.legend()
    plt.savefig(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--distance_tree", help="Input distance tree file", required=True)
    parser.add_argument("--simu_tree", help="Input distance file", required=True)
    parser.add_argument("--output", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.distance_tree, args.simu_tree, args.output)
