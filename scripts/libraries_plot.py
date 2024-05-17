import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.use('Agg')
hist_filled = {'alpha': 0.3, 'histtype': 'stepfilled'}
hist_step = {'histtype': 'step'}
my_dpi = 128
fontsize = 24
fontsize_legend = 18


def colors(x):
    cs = ["#D55E00", "#0072B2", "#F0E442", "#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#000000"]
    return [cs[idx % len(cs)] for idx in range(len(x))]


def hist_plot(x_input, x_label, output, xscale="log"):
    x = {k: np.array(v) for k, v in x_input.items()}
    if xscale == "log":
        xy_filtered = {k: (np.isfinite(x[k]) & (x[k] > 0)) for k in x.keys()}
    else:
        xy_filtered = {k: np.isfinite(x[k]) for k in x.keys()}
    x = {k: x[k][xy_filtered[k]] for k in x.keys() if len(x[k][xy_filtered[k]]) > 0}
    if len(x) == 0:
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
        ax.plot((x_mean, x_mean), (0, max_y), linewidth=3, color=color_models[id_m],
                label=f'{m.replace("_", " ").capitalize()} (mean {x_mean:.2g})')

    ax.set_xlabel(x_label, fontsize=fontsize)
    if xscale == "log":
        ax.set_xlim((0.95 * min_x, 1.05 * max_x))
        ax.set_xscale("log")
    elif xscale == "linear":
        ax.set_xlim((min_x, max_x))
        ax.set_xscale("linear")
    else:
        ax.set_xlim((-0.01, 1.01))
        ax.set_yscale("log")
    ax.set_ylim((0, max_y))
    ax.set_ylabel("Density", fontsize=fontsize)
    ax.legend(fontsize=fontsize_legend)
    plt.tight_layout()
    plt.savefig(output, format="pdf")
    plt.clf()
    print(output)
