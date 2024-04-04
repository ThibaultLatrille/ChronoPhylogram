import os
import argparse
from collections import defaultdict
import pandas as pd
from statsmodels.formula.api import ols
from libraries import *
import matplotlib.pyplot as plt


def main(distance_tree, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    simulator_dict = defaultdict(list)
    for file in distance_tree:
        simulator = os.path.basename(os.path.dirname(file))
        simulator_dict[simulator].append(open_tree(file, format_ete3=1))
    assert len(simulator_dict) == 2
    simu1, simu2 = list(simulator_dict.keys())
    assert len(simulator_dict[simu1]) == len(simulator_dict[simu2])
    n_seeds = len(simulator_dict[simu1])
    fig, axes = plt.subplots(n_seeds, n_seeds, figsize=(6 * n_seeds, 6 * n_seeds))
    for i, seed_i in enumerate(simulator_dict[simu1]):
        for j, seed_j in enumerate(simulator_dict[simu2]):
            x_list, y_list = [], []
            for node1, node2 in zip(seed_i.traverse(), seed_j.traverse()):
                assert node1.name == node2.name
                x_list.append(node1.dist)
                y_list.append(node2.dist)
            ax = axes[i, j]
            ax.scatter(x_list, y_list)
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
    plt.savefig(output)
    plt.close("all")
    plt.clf()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--distance_tree", help="Input distance tree file", nargs="+", required=True)
    parser.add_argument("--output", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.distance_tree, args.output)
