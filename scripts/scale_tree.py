import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from libraries import *


def parse_tree(tree_path: str) -> Tree:
    return rename_tree(prune_tree(open_tree(tree_path, format_ete3=1)))


def plot_tree(x_tree, y_tree, x_label, y_label, ax):
    distances_x, distances_y = [], []
    dict_x = {node.name: node.dist for node in x_tree.traverse()}
    dict_y = {node.name: node.dist for node in y_tree.traverse()}
    intersection = set(dict_x.keys()).intersection(dict_y.keys())
    for node_name in intersection:
        distances_x.append(dict_x[node_name])
        distances_y.append(dict_y[node_name])
    ax.scatter(distances_x, distances_y, alpha=0.4)
    # regression line and r-squared
    m, b = np.polyfit(distances_x, distances_y, 1)
    r2 = np.corrcoef(distances_x, distances_y)[0, 1] ** 2
    ax.plot(distances_x, m * np.array(distances_x) + b, color="red", label=f"y = {m:.2f}x + {b:.2f} (R² = {r2:.2f})")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title("Distance between the two trees")
    ax.legend()


def main(tree_1, tree_2, tree_output_1, tree_output_2):
    """
    Scale the tree_1 and tree_2 to have the same species and the same scale
    :return:
    """
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    x_label = os.path.basename(tree_2).replace(".tree", "")
    y_label = os.path.basename(tree_1).replace(".tree", "")

    t_1 = parse_tree(tree_1)
    t_2 = parse_tree(tree_2)
    intersection = None
    if set(t_1.get_leaf_names()) != set(t_2.get_leaf_names()):
        intersection = list(set(t_1.get_leaf_names()).intersection(set(t_2.get_leaf_names())))
        print(f"Intersection between tree_1 ({len(t_1)}) and tree_2 ({len(t_2)}) has {len(intersection)} leaves")
    t_1 = prune_tree(t_1, intersection)
    t_2 = prune_tree(t_2, intersection)
    plot_tree(t_1, t_2, x_label, y_label, axes[0])

    t_1 = scale_tree(t_1)
    t_2 = scale_tree(t_2)
    plot_tree(t_1, t_2, "Scaled " + x_label, "Scaled " + y_label, axes[1])
    assert len(t_1.get_leaf_names()) == len(t_2.get_leaf_names())
    t_1.write(outfile=tree_output_1, format=1)
    t_2.write(outfile=tree_output_2, format=1)
    # Plot the distance between the two trees as a scatter plot
    plt.savefig(tree_output_1 + ".pdf")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree_1", help="Input tree 1", required=True)
    parser.add_argument("--tree_2", help="Input tree 2", required=True)
    parser.add_argument("--tree_output_1", help="Output tree 1", required=True)
    parser.add_argument("--tree_output_2", help="Output tree 2", required=True)
    args = parser.parse_args()
    main(args.tree_1, args.tree_2, args.tree_output_1, args.tree_output_2)
