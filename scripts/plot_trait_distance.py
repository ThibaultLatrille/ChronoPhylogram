import argparse
import os.path
import numpy as np
from collections import defaultdict
from libraries import *
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

my_dpi = 96
fontsize = 24
fontsize_legend = 18


def saturation_fit(x, a, b):
    return a * (x / (b + x))


def linear_fit(x, a, b):
    return a * x + b


def polyfit(x_array, y_array, func, label, ax, color):
    popt, pcov = curve_fit(func, x_array, y_array, bounds=(0, np.inf))
    y_pred = func(x_array, *popt)
    r2 = r2_score(y_array, y_pred)
    if label == "saturation":
        reg = f"y = {popt[0]:.2g}x / ({popt[1]:.2g} + x) ($r^2$={r2:.2g})"
    else:
        reg = f"y = {popt[0]:.2g}x + {popt[1]:.2g} ($r^2$={r2:.2g})"
    x_line = np.linspace(np.min(x_array), np.max(x_array), 100)
    ax.plot(x_line, func(x_line, *popt), linestyle="--", label=reg, color=color)


def distance_plot(x_dico, y_dico, output, x_label='Tree distance', y_label='Phenotypic distance'):
    x_mean = {k: np.mean(v) for k, v in x_dico.items()}
    y_mean = {k: np.mean(v) for k, v in y_dico.items()}
    y_lower = {k: np.percentile(v, 10) for k, v in y_dico.items()}
    y_upper = {k: np.percentile(v, 90) for k, v in y_dico.items()}
    fig = plt.figure(figsize=(640 / my_dpi, 480 / my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(1, 1, 1)
    x_array = np.array(list(x_mean.values()))
    y_array = np.array(list(y_mean.values()))
    polyfit(x_array, y_array, saturation_fit, "saturation", ax, "blue")
    polyfit(x_array, y_array, linear_fit, "linear", ax, "red")
    ax.scatter(x_array, y_array, s=12, color="black")
    for k in x_mean.keys():
        ax.errorbar(x_mean[k], y_mean[k], color="black", alpha=0.05, linewidth=0.5,
                    yerr=[[y_mean[k] - y_lower[k]], [y_upper[k] - y_mean[k]]])
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    ax.legend(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.close('all')
    print(output)


def scatter_plot(data, ax: plt.Axes, title: str, x_label: str, y_label: str, average: bool = False):
    r2_list = []
    for key, (x, y, y_std) in data.items():
        if y_std is not None:
            # Display the range of the inference nodes, hide the center points
            ax.errorbar(x, y, yerr=y_std, fmt='none', alpha=0.5, zorder=-1)
        assert len(x) == len(y), "Different number of x and y values"
        # regression line and r-squared
        a, b = np.polyfit(x, y, 1)
        r2 = np.corrcoef(x, y)[0, 1] ** 2
        x_lim = np.array([min(0, min(x)), max(x)])
        if average:
            ax.scatter(x, y, alpha=0.25)
        else:
            label = f"{key} (R² = {r2:.2f}, n = {len(x)})"
            ax.scatter(x, y, label=label)
        ax.plot(x_lim, a * x_lim + b, color=ax.collections[-1].get_facecolor()[0], linestyle="--")
        r2_list.append(r2)
    ax.axhline(0, color="black", linestyle="--")
    ax.set_xlabel(x_label, fontsize=fontsize)
    ax.set_ylabel(y_label, fontsize=fontsize)
    if average:
        title += f" (mean R² = {np.mean(r2_list):.2f})"
    else:
        ax.legend(fontsize=fontsize_legend)
    ax.set_title(title, fontsize=fontsize)


def main(input_distance_tree: str, input_annotated_tree_list: list[str], output_path: str):
    d_tree = open_tree(input_distance_tree, format_ete3=1)
    dict_d = {node.name: node for node in d_tree.traverse()}
    dict_x, dict_y = defaultdict(list), defaultdict(list)
    for input_annotated_tree in input_annotated_tree_list:
        annotated_tree = open_tree(input_annotated_tree, format_ete3=1)
        dict_annotated = {node.name: node for node in annotated_tree.traverse()}
        intersection = set(dict_d.keys()).intersection(dict_annotated.keys())
        for node_name in intersection:
            n_d = dict_d[node_name]
            n_t = dict_annotated[node_name]
            if n_t.is_root() or n_t.up.is_root():
                continue
            assert n_t.name == n_d.name, f"Taxon names do not match: {n_t.name} != {n_d.name}"
            dict_x[n_t.name].append(np.sqrt(n_d.dist))
            trait = float(getattr(n_t, "Phenotype_mean"))
            trait_parent = float(getattr(n_t.up, "Phenotype_mean"))
            dict_y[n_t.name].append(np.abs(trait - trait_parent))

    leaves = set([node for node in dict_x.keys() if dict_d[node].is_leaf()])
    nodes = set([node for node in dict_x.keys() if not dict_d[node].is_leaf()])

    x_label = "$\\sqrt{" + ("d" if "Phylo" in os.path.basename(input_distance_tree) else "T") + "}$"
    y_label = "|Δz|"
    fig, ax = plt.subplots(1, 1, figsize=(640 / my_dpi, 480 / my_dpi), dpi=my_dpi)
    x_mean = {k: np.mean(v) for k, v in dict_x.items()}
    y_mean = {k: np.mean(v) for k, v in dict_y.items()}
    y_std = {k: np.std(v) for k, v in dict_y.items()}
    dico_data = {"all": (list(x_mean.values()), list(y_mean.values()), list(y_std.values())),
                 "leaves": ([x_mean[k] for k in leaves], [y_mean[k] for k in leaves], [y_std[k] for k in leaves]),
                 "nodes": ([x_mean[k] for k in nodes], [y_mean[k] for k in nodes], [y_std[k] for k in nodes])}
    scatter_plot(dico_data, ax, "Average across simulations", x_label=x_label, y_label=y_label)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close('all')

    # Saturation plot
    distance_output = replace_last(output_path, ".pdf", "_saturation.pdf")
    distance_plot(dict_x, dict_y, distance_output, x_label=x_label, y_label=y_label)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--distance_tree", help="Input distance tree file", required=True)
    parser.add_argument("--annotated_tree", help="Input annotated tree file", nargs="+", required=True)
    parser.add_argument("--output", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.distance_tree, args.annotated_tree, args.output)
