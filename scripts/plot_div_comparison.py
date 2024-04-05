import os
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
from statsmodels.formula.api import ols
from libraries import *
import matplotlib.pyplot as plt


def tree_to_dict(tree):
    tree_dict = {}
    for node in tree.traverse():
        tree_dict[node.name] = node.dist
    return tree_dict


def main(distance_tree, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    simulator_dict = defaultdict(list)
    for file in distance_tree:
        simulator = os.path.basename(os.path.dirname(file))
        simulator_dict[simulator].append(tree_to_dict(open_tree(file, format_ete3=1)))
    simu_list = list(simulator_dict.keys())
    simu1 = simu_list[0]
    for id_simu, simu2 in enumerate(simu_list[1:]):
        assert len(simulator_dict[simu1]) == len(simulator_dict[simu2])
        n_seeds = min(len(simulator_dict[simu1]), 10)
        list_simu1 = np.random.choice(simulator_dict[simu1], n_seeds, replace=False)
        list_simu2 = np.random.choice(simulator_dict[simu2], n_seeds, replace=False)
        fig, axes = plt.subplots(n_seeds, n_seeds, figsize=(6 * n_seeds, 6 * n_seeds))
        for i, dict_i in enumerate(list_simu1):
            for j, dict_j in enumerate(list_simu2):
                x_list, y_list = [], []
                intersection = set(dict_i.keys()).intersection(dict_j.keys())
                for node_name in intersection:
                    x_list.append(dict_i[node_name])
                    y_list.append(dict_j[node_name])
                ax = axes[i, j]
                ax.scatter(x_list, y_list, alpha=0.4)
                ols_res = ols(f"y ~ x", data=pd.DataFrame({"x": x_list, "y": y_list})).fit()
                title = f"slope: {ols_res.params['x']:.2f}, r2: {ols_res.rsquared:.2f}"
                # plot regression line
                x_min, x_max = min(x_list), max(x_list)
                ax.plot([x_min, x_max], [ols_res.params["x"] * x_min + ols_res.params["Intercept"],
                                         ols_res.params["x"] * x_max + ols_res.params["Intercept"]], color="red")
                # plot y=x line
                ax.plot([x_min, x_max], [x_min, x_max], color="black", linestyle="--")
                ax.set_title(title)
                ax.set_xlabel(simu1)
                ax.set_ylabel(simu2)
        fig.tight_layout()
        if id_simu == 0:
            plt.savefig(output)
        else:
            plt.savefig(output.replace(".pdf", f"_{id_simu}.pdf"))
        plt.close("all")
        plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--distance_tree", help="Input distance tree file", nargs="+", required=True)
    parser.add_argument("--output", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.distance_tree, args.output)
