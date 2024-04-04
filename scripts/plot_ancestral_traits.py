import os
import argparse
import numpy as np
from libraries import *
import matplotlib.pyplot as plt


def plot_ancestral_traits(input_inference_tree: str, input_simu_tree: str, ax):
    inference_tree = open_tree(input_inference_tree, format_ete3=1)
    simu_tree = open_tree(input_simu_tree, format_ete3=1)
    simu_nodes, inference_nodes_min, inference_nodes, inference_nodes_max = [], [], [], []
    for simu_node, inference_node in zip(simu_tree.traverse(), inference_tree.traverse()):
        if simu_node.is_root() or inference_node.is_leaf():
            continue
        assert simu_node.name == inference_node.name, f"Taxon names do not match: {simu_node.name} != {inference_node.name}"
        simu_nodes.append(float(getattr(simu_node, "Phenotype_mean")))

        inference_nodes_min.append(float(getattr(inference_node, "Phenotype_mean_min")))
        inference_nodes.append(float(getattr(inference_node, "Phenotype_mean")))
        inference_nodes_max.append(float(getattr(inference_node, "Phenotype_mean_max")))
    inference_nodes = np.array(inference_nodes)
    inference_nodes_min = np.array(inference_nodes_min)
    inference_nodes_max = np.array(inference_nodes_max)
    assert len(simu_nodes) == len(inference_nodes), "Different number of x and y values"
    # Scatter plot
    ax.scatter(simu_nodes, inference_nodes)
    # Display the range of the inference nodes
    ax.errorbar(simu_nodes, inference_nodes,
                 yerr=[inference_nodes - inference_nodes_min, inference_nodes_max - inference_nodes],
                 fmt='o', color='black', alpha=0.5)
    # regression line and r-squared
    m, b = np.polyfit(simu_nodes, inference_nodes, 1)
    r2 = np.corrcoef(simu_nodes, inference_nodes)[0, 1] ** 2
    ax.plot(simu_nodes, m * np.array(simu_nodes) + b, color="red", label=f"y = {m:.2f}x + {b:.2f} (R² = {r2:.2f})")
    # Plot the identity line
    ax.plot([min(simu_nodes), max(simu_nodes)], [min(simu_nodes), max(simu_nodes)], color="black", linestyle="--",
             label="y=x")
    ax.set_xlabel("Phenotype_mean (simu)")
    ax.set_ylabel("Phenotype_mean (inference)")
    ax.set_title(os.path.basename(input_simu_tree.replace(".nhx.gz", "")))
    ax.legend()
    return r2


def main(input_inference_tree_list: str, input_simu_tree_list: str, output_path: str):
    input_simu_tree_list = sorted(input_simu_tree_list)
    input_inference_tree_list = sorted(input_inference_tree_list)
    assert len(input_simu_tree_list) == len(input_inference_tree_list), "Different number of input files"
    ncols = int(np.sqrt(len(input_simu_tree_list)) + 1)
    nrows = int(len(input_simu_tree_list) / ncols + 1)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 6, nrows * 6))
    i = 1
    r2_list = []
    for input_inference_tree, input_simu_tree in zip(input_inference_tree_list, input_simu_tree_list):
        seed_simu = int(os.path.basename(input_simu_tree).split("_")[1].replace(".nhx.gz", "").replace("seed", ""))
        seed_inference = int(os.path.basename(input_inference_tree).split("_")[2].replace(".Phenotype", "").replace("seed", ""))
        assert seed_simu == seed_inference, "Different seeds"
        ax = axes[i // ncols, i % ncols]
        r2 = plot_ancestral_traits(input_inference_tree, input_simu_tree, ax)
        r2_list.append(r2)
        i += 1
    ax = axes[0, 0]
    # plot the histogram of R² values
    ax.hist(r2_list, bins=20)
    ax.set_xlabel("R²")
    ax.set_ylabel("Frequency")
    ax.set_title("Histogram of R² values")
    ax.set_xlim(0, 1)
    ax.vlines(np.mean(r2_list), 0, max(ax.get_ylim()), color="red", label=f"Mean R² = {np.mean(r2_list):.2f}")
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--inference_tree", help="Input inference tree file", nargs="+", required=True)
    parser.add_argument("--simu_tree", help="Input simulation file", nargs="+", required=True)
    parser.add_argument("--output", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.inference_tree, args.simu_tree, args.output)
