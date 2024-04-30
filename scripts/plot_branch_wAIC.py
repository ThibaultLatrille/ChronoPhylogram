import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from libraries import *
from merge_wAIC_simulations import parse_bayescode_list


def tree_to_dict(tree_path: str) -> dict:
    t = open_tree(tree_path, format_ete3=1)
    return {node.name: node.dist for node in t.traverse()}


def plot_tree(dict_x, dict_y, dict_z, x_label, y_label, z_label, title, ax):
    distances_x, distances_y, z_values = [], [], []
    intersection = set(dict_x.keys()).intersection(dict_y.keys()).intersection(dict_z.keys())
    for node_name in intersection:
        distances_x.append(dict_x[node_name])
        distances_y.append(dict_y[node_name])
        z_values.append(dict_z[node_name])
    # Color by z_values, where z_values of 0 are displayed as white, and the colormap centered on 0
    cmap = plt.colormaps.get_cmap("coolwarm")
    vmin, vmax = min(z_values), max(z_values)
    # Center the colormap on 0
    if abs(vmin) > abs(vmax) and vmin < 0 < vmax:
        vmax = abs(vmin)
    else:
        vmin = -abs(vmax)
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    sc = ax.scatter(distances_x, distances_y, c=z_values, cmap=cmap, norm=norm, label=f"n={len(z_values)}")
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(z_label)

    # Plot identity line
    ax.plot([min(distances_x), max(distances_x)], [min(distances_x), max(distances_x)], color="black", linestyle="--")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend()


def main(tree_x, tree_y, bayescode_list, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    x_label = os.path.basename(tree_x).replace(".tree", "")
    y_label = os.path.basename(tree_y).replace(".tree", "")
    dict_x = tree_to_dict(tree_x)
    dict_y = tree_to_dict(tree_y)

    df = parse_bayescode_list(bayescode_list, branch=True)
    n = len(df["experiment"].unique())
    ncols = int(np.sqrt(n) + 1)
    nrows = int(n / ncols + 1)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 8, nrows * 5))
    average_waic = df.groupby("name")["ΔwAIC"].mean()
    dico_waic = {k: v for k, v in zip(average_waic.index, average_waic)}
    ax = axes[0, 0]
    plot_tree(dict_x, dict_y, dico_waic, x_label, y_label, "ΔwAIC", "Average", ax)
    i = 1
    for experiment, df_experiment in df.groupby("experiment"):
        dico_waic = {}
        for row_index, row in df_experiment.iterrows():
            dico_waic[row["name"]] = row["ΔwAIC"]
        plot_tree(dict_x, dict_y, dico_waic, x_label, y_label, "ΔwAIC", experiment, axes[i // ncols, i % ncols])
        i += 1
    plt.tight_layout()
    plt.savefig(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree_x", help="Input tree x", required=True)
    parser.add_argument("--tree_y", help="Input tree y", required=True)
    parser.add_argument("--bayescode", nargs="+", help="Input BayesCode file", required=True)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.tree_x, args.tree_y, args.bayescode, args.output)
