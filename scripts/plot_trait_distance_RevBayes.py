import argparse
import os.path
import numpy as np
from collections import defaultdict
from libraries import *
import matplotlib.pyplot as plt
from libraries_plot import cs_simu_models, gr_simu_models

my_dpi = 96
fontsize = 24
fontsize_legend = 18


def scatter_plot(data, ax: plt.Axes, title: str, x_label: str, y_label: str, average: bool = False, color: str = "black"):
    r2_list = []
    for key, (x, y, y_std) in data.items():
        if y_std is not None:
            # Display the range of the inference nodes, hide the center points
            ax.errorbar(x, y, yerr=y_std, fmt='none', alpha=0.1, zorder=-1, ecolor=color)
        assert len(x) == len(y), "Different number of x and y values"
        # regression line and r-squared
        a, b = np.polyfit(x, y, 1)
        r2 = np.corrcoef(x, y)[0, 1] ** 2
        x_lim = np.array([min(0, min(x)), max(x)])
        if average:
            ax.scatter(x, y, alpha=0.25, color=color)
        else:
            ax.scatter(x, y, label=f"{key} (n = {len(x)})", color=color)
        ax.plot(x_lim, a * x_lim + b, color="black", linestyle="--", label=f"Linear regression (R² = {r2:.2f})")
        r2_list.append(r2)
    # ax.axhline(0, color="black", linestyle="--")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    if average:
        title += f" (mean R² = {np.mean(r2_list):.2f})"
    else:
        ax.legend(fontsize=fontsize_legend)
    ax.set_title(title, fontsize=fontsize)


def main(input_distance_tree: str, input_annotated_tree_list: list[str], output_path: str):
    d_tree = open_tree(input_distance_tree, format_ete3=1)
    d_tree.name = "Root"
    dict_d = {node.name: node for node in d_tree.traverse()}
    dict_x, dict_y = defaultdict(list), defaultdict(list)
    for input_annotated_tree in input_annotated_tree_list:
        annotated_tree = open_tree(input_annotated_tree, format_ete3=1)
        dict_annotated = {node.name: node for node in annotated_tree.traverse()}
        for d_node in d_tree.traverse():
            if d_node.name not in dict_annotated.keys():
                # find the same node in the inference tree by matching the leaves names of the subtree
                leaves_names = set(d_node.get_leaf_names())
                for annot_node in annotated_tree.traverse():
                    if set(annot_node.get_leaf_names()) == leaves_names:
                        dict_annotated[d_node.name] = annot_node
                        annot_node.name = d_node.name
                        break
        intersection = set(dict_d.keys()).intersection(dict_annotated.keys())
        for node_name in intersection:
            n_d = dict_d[node_name]
            n_t = dict_annotated[node_name]
            if n_t.is_root() or n_t.up.is_root() or n_d.dist == 0:
                continue
            assert n_t.name == n_d.name, f"Taxon names do not match: {n_t.name} != {n_d.name}"
            dict_x[n_t.name].append(np.sqrt(n_d.dist))
            trait = float(getattr(n_t, "Phenotype_mean"))
            trait_parent = float(getattr(n_t.up, "Phenotype_mean"))
            dict_y[n_t.name].append(np.abs(trait - trait_parent))

    x_label = "$\\sqrt{" + ("d" if "Phylo" in os.path.basename(input_distance_tree) else "T") + "}$"
    y_label = "|Δz|"
    fig, ax = plt.subplots(1, 1, figsize=(640 / my_dpi, 480 / my_dpi), dpi=my_dpi)
    x_mean = {k: np.mean(v) for k, v in dict_x.items()}
    y_mean = {k: np.mean(v) for k, v in dict_y.items()}
    y_std = {k: np.std(v) for k, v in dict_y.items()}
    dico_data = {"All branches": (list(x_mean.values()), list(y_mean.values()), list(y_std.values()))}
    color = cs_simu_models[gr_simu_models(os.path.basename(output_path))][0]
    scatter_plot(dico_data, ax, "", x_label=x_label, y_label=y_label, color=color)
    plt.tight_layout()
    plt.savefig(output_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--distance_tree", help="Input distance tree file", required=True)
    parser.add_argument("--annotated_tree", help="Input annotated tree file", nargs="+", required=True)
    parser.add_argument("--output", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.distance_tree, args.annotated_tree, args.output)
