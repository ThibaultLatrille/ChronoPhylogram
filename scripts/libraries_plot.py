import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pyparsing import White

mpl.use('Agg')
hist_filled = {'alpha': 0.3, 'histtype': 'stepfilled'}
hist_step = {'histtype': 'step'}
my_dpi = 128
fontsize = 24
fontsize_legend = 18
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
    cs = {"N": (ORANGE, "Neutral"), "B": (BLUE, "Multiple optimal"), "S": (RED, "Moving optimum")}
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
    # Put the neutral model first, then the multiple optimal model, then the moving optimum model
    return list(sorted(list_x, key=lambda x: 0 if "neutral" in x else (1 if "multi" in x and "optimum" in x else 2)))


def vert_boxplot(x_input, y_label, output, yscale="linear", format_label=None, empirical=False):
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
    data = [x[m] for m in sorted_x]
    sns.boxplot(data=data, ax=ax, palette=color_models, fliersize=0, log_scale=(yscale == "log"))
    # Display points on top of the boxplot
    sns.swarmplot(data=data, ax=ax, color="black", size=1, edgecolor="black", linewidth=1)
    ax.set_ylabel(y_label, fontsize=fontsize_legend)
    labels = [m.replace("_", " ") for m in sorted_x]
    if format_label is not None:
        labels = [format_label(l) for l in labels]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=fontsize_legend)
    ax.legend(handles=handles, fontsize=fontsize_legend)
    plt.xticks(rotation=45, ha='right')
    if yscale == "uniform":
        ax.set_ylim((-0.01, 1.01))
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.clf()
    plt.close("all")
    print(output)
