import numpy as np
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

mpl.use('Agg')
hist_filled = {'alpha': 0.3, 'histtype': 'stepfilled'}
hist_step = {'histtype': 'step'}
my_dpi = 128
fontsize = 24
fontsize_legend = 18


def colors(x):
    cs = ["#D55E00", "#0072B2", "#F0E442", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#000000"]
    return [cs[idx % len(cs)] for idx in range(len(x))]


def color_simu_models(x):
    cs = {"N": "#F0E442", "S": "#E69F00"}
    color_list = [(cs["N"] if "neutral" in xl.lower() else cs["S"]) for xl in x.keys()]
    handles = [Rectangle((0, 0), 1, 1, color=cs["N"], ec="k", lw=1, label="Neutral"),
               Rectangle((0, 0), 1, 1, color=cs["S"], ec="k", lw=1, label="Moving optimum")]
    return color_list, handles


def filter_x(x_input, xscale):
    x = {k: np.array(v) for k, v in x_input.items()}
    if xscale == "log":
        xy_filtered = {k: (np.isfinite(x[k]) & (x[k] > 0)) for k in x.keys()}
    else:
        xy_filtered = {k: np.isfinite(x[k]) for k in x.keys()}
    return {k: x[k][xy_filtered[k]] for k in x.keys() if len(x[k][xy_filtered[k]]) > 0}


def hist_plot(x_input, x_label, output, xscale="log", format_label=None):
    x = filter_x(x_input, xscale)
    if len(x) <= 1:
        return
    color_models = colors(x)
    fig = plt.figure(figsize=(1280 / my_dpi, 640 / my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(1, 1, 1)
    min_x = min([v.min() for v in x.values()])
    max_x = max([v.max() for v in x.values()])
    if xscale == "log":
        min_x = max(1e-6, min_x)
        bins = np.geomspace(min_x, max_x, 100)
    else:
        bins = np.linspace(min_x, max_x, 50)
    hist, _, _ = ax.hist(x.values(), bins=bins, color=color_models, **hist_filled)
    hist, _, _ = ax.hist(x.values(), bins=bins, color=color_models, **hist_step)
    max_y = 1.2 * (np.max([np.max(h) for h in hist]) if len(x) > 1 else np.max(hist))

    for id_m, m in enumerate(x):
        x_mean = np.mean(x[m])
        label = m.replace("_", " ").capitalize()
        if format_label is not None:
            label = format_label(label)
        ax.plot((x_mean, x_mean), (0, max_y), linewidth=3, color=color_models[id_m], label=f"{label}: {x_mean:.2f}")

    ax.set_xlabel(x_label, fontsize=fontsize)
    if xscale == "log":
        ax.set_xlim((0.95 * min_x, 1.05 * max_x))
        ax.set_xscale("log")
        ax.set_ylim((0, max_y))
    elif xscale == "linear":
        ax.set_xlim((min_x, max_x))
        ax.set_xscale("linear")
        ax.set_ylim((0, max_y))
    else:
        ax.set_xlim((-0.01, 1.01))
        ax.set_yscale("log")
    ax.set_ylabel("Density", fontsize=fontsize)
    ax.legend(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.clf()
    plt.close("all")
    print(output)


def vert_boxplot(x_input, x_label, output, xscale="linear", format_label=None):
    x = filter_x(x_input, xscale)
    if len(x) <= 1:
        return
    color_models, handles = color_simu_models(x)
    fig = plt.figure(figsize=(1280 / my_dpi, 640 / my_dpi), dpi=my_dpi)
    ax = fig.add_subplot(1, 1, 1)
    sns.boxplot(data=[x[m] for m in x], ax=ax, palette=color_models, fliersize=0, log_scale=(xscale == "log"))
    # Display points on top of the boxplot
    sns.swarmplot(data=[x[m] for m in x], ax=ax, color="black")
    ax.set_ylabel(x_label, fontsize=fontsize_legend)
    labels = [m.replace("_", " ") for m in x]
    if format_label is not None:
        labels = [format_label(l) for l in labels]
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=fontsize_legend)
    ax.legend(handles=handles, fontsize=fontsize_legend)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.clf()
    plt.close("all")
    print(output)
