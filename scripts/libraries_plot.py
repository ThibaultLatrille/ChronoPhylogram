import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from statannotations.Annotator import Annotator

mpl.use('Agg')
hist_filled = {'alpha': 0.3, 'histtype': 'stepfilled'}
hist_step = {'histtype': 'step'}
my_dpi = 128
fontsize = 24
fontsize_legend = 16
RED = "#D55E00"
BLUE = "#0072B2"
GREEN = "#009E73"
ORANGE = "#E69F00"
PURPLE = "#CC79A7"
YELLOW = "#F0E442"
CYAN = "#56B4E9"
WHITE = "#FFFFFF"
BLACK = "#000000"


def colors(x):
    cs = [RED, BLUE, YELLOW, ORANGE, CYAN, GREEN, PURPLE, BLACK]
    return [cs[idx % len(cs)] for idx in range(len(x))]


def gr_simu_models(x):
    return "N" if ("neutral" in x) else ("B" if (("multi" in x) and ("optimum" in x)) else "S")


def color_simu_models(x_list):
    cs = {"N": (ORANGE, "Neutral"), "B": (BLUE, "Multiple optima"), "S": (RED, "Moving optimum")}
    grs = [gr_simu_models(xl) for xl in x_list]
    color_list = [cs[g][0] for g in grs]
    handles = [Rectangle((0, 0), 1, 1, color=v[0], ec="k", lw=1, label=v[1]) for c, v in cs.items() if c in set(grs)]
    return color_list, handles


def color_mf_models(x_list):
    cs = {"m": (RED, "Males"), "f": (BLUE, "Females"), "b": (ORANGE, "Both")}
    grs = [xl[0] for xl in x_list]
    color_list = [cs[g][0] for g in grs]
    handles = [Rectangle((0, 0), 1, 1, color=v[0], ec="k", lw=1, label=v[1]) for c, v in cs.items() if c in set(grs)]
    return color_list, handles


def color_empirical_models(x_list):
    cs = {"bodyMass": (BLUE, "Body Mass"), "brainMass": (GREEN, "Brain mass")}
    grs = [("bodyMass" if "body" in xl.lower() else "brainMass") for xl in x_list]
    color_list = [cs[g][0] for g in grs]
    handles = [Rectangle((0, 0), 1, 1, color=v[0], ec="k", lw=1, label=v[1]) for c, v in cs.items() if c in set(grs)]
    return color_list, handles


def filter_x(x_input, xscale):
    x = {k: np.array(v) for k, v in x_input.items()}
    if xscale == "log":
        xy_filtered = {k: (np.isfinite(x[k]) & (x[k] > 0)) for k in x.keys()}
    else:
        xy_filtered = {k: np.isfinite(x[k]) for k in x.keys()}
    return {k: x[k][xy_filtered[k]] for k in x.keys() if len(x[k][xy_filtered[k]]) > 0}


def sort_x(list_x):
    # Put the neutral model first, then the multiple optima model, then the moving optimum model
    return list(sorted(list_x, key=lambda x: 0 if "neutral" in x else (1 if "multi" in x and "optimum" in x else 2)))


def remove_chrono_phylo(x):
    return x.replace("Chrono", "").replace("Phylo", "")


def vert_boxplot(x_input, y_label, output, yscale="linear", format_label=None, empirical=False, prior=None):
    x = filter_x(x_input, yscale)
    if len(x) < 1:
        return
    sorted_x = sort_x(x.keys())
    if empirical:
        color_models, handles = color_empirical_models(sorted_x)
    else:
        color_models, handles = color_simu_models(sorted_x)
    fig = plt.figure(figsize=(1280 / my_dpi, 640 / my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(1, 1, 1)
    df = pd.DataFrame({"y": np.concatenate([x[k] for k in sorted_x]),
                       "x": np.concatenate([[k] * len(x[k]) for k in sorted_x])})
    sns.violinplot(data=df, x="x", y="y", ax=ax, palette=color_models, log_scale=(yscale == "log"), inner="stick",
                   cut=0, legend=False, hue="x")
    ax.set_xlabel("")
    if "".join(sorted_x).count("Both") > 0:
        handles = []
        for i, k in enumerate(sorted_x):
            # Print the proportion above the 95% and below the 5% quantiles
            count_95 = len([1 for v in x[k] if v > .95])
            count_05 = len([1 for v in x[k] if v < .05])
            if count_95 > 0:
                ax.text(i, 0.95, f"{count_95} / {len(x[k])} > 0.95",
                        bbox=dict(facecolor="white", edgecolor="black", boxstyle="round", pad=0.15),
                        fontsize=13, ha="center", va="center", zorder=10)
            if count_05 > 0:
                ax.text(i, 0.05, f"{count_05} / {len(x[k])} < 0.05",
                        bbox=dict(facecolor="white", edgecolor="black", boxstyle="round", pad=0.15),
                        fontsize=13, ha="center", va="center", zorder=10)
            # h.set_label(r"$\hat{\mu}=$" + f"{np.mean(x[k]):.2f}")
            ax.text(i, 0.55, r"$\hat{\mu} = $" + str(round(np.mean(x[k]), 2)),
                    bbox=dict(facecolor="white", edgecolor=color_models[i], boxstyle="round", pad=0.4),
                    fontsize=13, ha="center", va="center", zorder=10)
        if yscale == "uniform":
            ax.axhline(0.05, color="grey", linestyle="--", linewidth=1)
            ax.axhline(0.95, color="grey", linestyle="--", linewidth=1)
    else:
        # Add the mean above each violin
        for i, k in enumerate(sorted_x):
            x_mean = np.mean(x[k])
            ax.text(i, x_mean, r"$\hat{\mu} = $" + str(round(x_mean, 3 if x_mean < 0.01 else 2)),
                    bbox=dict(facecolor="white", edgecolor="black", boxstyle="round", pad=0.15),
                    fontsize=13, ha="center", va="center", zorder=10)
        # Test the difference between the models
        pairs = []
        # They must be different only for chronogram/phylogram in the name
        for i in range(1, len(sorted_x)):
            if remove_chrono_phylo(sorted_x[i - 1]) == remove_chrono_phylo(sorted_x[i]):
                pairs.append((sorted_x[i - 1], sorted_x[i]))
        if len(pairs) > 0:
            annotator = Annotator(ax, pairs, data=df, x="x", y="y", order=sorted_x)
            annotator.configure(test='Wilcoxon', text_format='simple', loc='outside', verbose=0)
            annotator.apply_and_annotate()
    ax.set_ylabel(y_label, fontsize=fontsize_legend)
    labels = [m.replace("_", " ") for m in sorted_x]
    if format_label is not None:
        labels = [format_label(l) for l in labels]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=fontsize_legend)
    if len(x) < 4:
        plt.xticks(rotation=0, ha="center")
    else:
        plt.xticks(rotation=45, ha="right")
    if prior is not None:
        ax.axhline(prior, color="black", linestyle="-", linewidth=1)
        handles.append(Line2D([0], [0], color="black", linestyle="-", label=f"Prior ({prior})"))
    if yscale == "uniform":
        ax.set_ylim((0.0, 1.0))
    if len(handles) > 1:
        ax.legend(handles=handles, fontsize=12)
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.clf()
    plt.close("all")
    print(output)
