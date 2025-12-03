import os
import argparse
from glob import glob
import numpy as np
from natsort import natsorted
import pandas as pd
from collections import defaultdict
from os.path import basename, isdir
from libraries_plot import vert_boxplot
import seaborn as sns
import matplotlib.pyplot as plt

my_dpi = 150
fontsize = 24
fontsize_legend = 18


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def plot_violin(df_out: pd.DataFrame, col: str, output: str):
    fig, axes = plt.subplots(2, 1, figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    ax = axes[0]
    df_out.sort_values(by=["seed", "gram"], inplace=True)
    palette = {"Chrono": sns.color_palette("tab10")[0], "Phylo": sns.color_palette("tab10")[1]}
    sns.violinplot(data=df_out, x="seed", y=col, hue="gram", split=True, inner="quart", fill=True,
                   density_norm='width', palette=palette, ax=ax)
    ax.legend(fontsize=fontsize_legend)
    ax.set_xticks(range(len(ax.get_xticklabels())))
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    ax = axes[1]
    for dataset, df_seed in df_out.groupby("seed"):
        phylo = df_seed[df_seed["gram"] == "Phylo"][col].values
        chrono = df_seed[df_seed["gram"] == "Chrono"][col].values
        assert len(phylo) == len(chrono), f"Length of phylo ({len(phylo)}) and chrono ({len(chrono)}) are different"
        diff = phylo - chrono
        sns.kdeplot(diff, label=dataset, ax=ax)
    ax.set_xlabel("lnprob(Phylo) - lnprob(Chrono)", fontsize=fontsize)
    ax.axvline(0, color='black', linestyle='--')
    plt.tight_layout()
    plt.savefig(output)
    plt.clf()
    plt.close("all")


def open_log(path):
    df = pd.read_csv(path, sep="\t")
    if "alpha" in df.columns:
        df["t_half"] = np.log(2) / df["alpha"]
    if "var_rates" in df.columns and "sigma2_root" in df.columns:
        df["std_rates"] = np.sqrt(df["var_rates"] / df["sigma2_root"])
    return df


def posterior(df, col):
    return np.nanmean(df[col])


def get_dicts(simu_models: list, replicates: dict, rb_models: list, parameters: dict):
    post_dict = defaultdict(lambda: defaultdict(list))
    output_rows = list()
    list_df = []

    for simu_m in simu_models:
        for folderpath in replicates[simu_m]:
            n, gram, seed = basename(folderpath).split("_")
            for rb_model in rb_models:
                trace_path = f"{folderpath}/{rb_model}.log.gz"
                if not os.path.exists(trace_path):
                    continue

                trace_df = open_log(trace_path)
                if rb_model in ["simple_BM_REML", "simple_BM_nodes"]:
                    df = pd.DataFrame({"Posterior": trace_df["Posterior"]})
                    df["model"] = rb_model
                    df["gram"] = gram
                    df["seed"] = seed
                    df["simu"] = simu_m
                    df["dataset"] = gram + "_" + seed
                    list_df.append(df)

                key_name = f"{simu_m}_{gram}_{rb_model.replace("simple_", "S").replace("relaxed_", "R")}"



                row = {"simu": simu_m, "gram": gram, "seed": seed, "model": rb_model}
                for col in trace_df.columns:
                    if col not in parameters:
                        continue
                    if "multi" in simu_m and (
                            col in ["is_OU", "var_multiplier", "num_rate_changes", "num_theta_changes"]):
                        continue
                    post_v = posterior(trace_df, col=col)
                    post_dict[col][key_name].append(post_v)
                    row[col] = post_v
                output_rows.append(row)
    # output_rows is a list of dicts (missing values for some parameters)
    df_post = pd.DataFrame(output_rows)
    return pd.concat(list_df), post_dict, df_post


def plot_trace(df_out: pd.DataFrame, col: str, output: str):
    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    for _, df in df_out.groupby(["seed", "gram"]):
        plt.plot(df[col])
    plt.savefig(output)
    plt.clf()


def format_label(l):
    rm_list = {"Both", "SBM", "RBM", "Switch", "RJ", "SOU", "ROU", "nodes", "REML"}
    for rm in rm_list:
        l = l.replace(rm, "")
    # Remove "neutral", "moving", "multi", "optimum" only there is no other word
    s = l.replace("neutral", "").replace("moving", "").replace("multi", "").replace("optimum", "")
    if len(s.strip()) == 0:
        s = l
        s = s.replace("multi_optimum", "multiple_optima")
    while "  " in s:
        s = s.replace("  ", " ")
    return s.replace("Chrono", "Chronogram").replace("Phylo", "Phylogram").title()


def main(folder, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    rb_models = ["relaxed_BM_RJ", "relaxed_OU_RJ", "simple_OU_RJ", "simple_BM_REML", "simple_BM_nodes",
                 "simple_BM_SwitchREML", "simple_BM_Switchnodes"]
    simu_model_prefs = {"moving_optimum": 0, "multi_optimum": 1, "neutral": 3}
    simu_models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p) and basename(p) in simu_model_prefs}
    simu_models = list(sorted(simu_models_path, key=lambda x: simu_model_prefs[x] if x in simu_model_prefs else -1))

    replicates = {m: natsorted([i for i in glob(f"{p}/*") if isdir(i)]) for m, p in simu_models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    parameters = {"num_rate_changes": ("Number of rate changes", "linear", 1.0, "n"),
                  "std_rates": ("std rates", "log", None, "n"),
                  "var_multiplier": ("Variance in the rates", "log", None, "v"),
                  "var_rates": ("Variance in the rates", "log", None, "v"),
                  "num_theta_changes": ("Number of optimum changes", "linear", 1.0, "m"),
                  "theta_multiplier": ("theta", "log", None, "n"),
                  "is_BM": ("Support for BM over OU", "uniform", 0.5, "p_{BM}"),
                  "is_OU": ("Support for OU over BM", "uniform", 0.5, "p_{OU}"),
                  "is_nuc": ("Support for a phylogram", "uniform", 0.5, "\\pi"),
                  "sigma": ("sigma", "log", None, "\\sigma"), "theta": ("theta", "linear", None, "\\theta"),
                  "alpha": ("alpha", "log", None, "\\alpha"), "t_half": ("t 1/2", "log", None, "T_{1/2}")}
    print(f"Importing data from {folder}")
    trace_df, post_dict, output_df = get_dicts(simu_models, replicates, rb_models, parameters)

    rename = lambda x: output.replace(".tsv", x)
    for col, dict_input in post_dict.items():
        y_label, yscale, prior, var_name = parameters[col]
        y_label = y_label + f" (${var_name}$)"
        if "_" in var_name:
            var_name_mean = f"\\overline{{{var_name.split('_')[0]}}}_{{{var_name.split('_')[1]}}}"
        else:
            var_name_mean = f"\\overline{{{var_name}}}"
        print(f"Plotting {col}")
        vert_boxplot(dict_input, y_label, rename(f".boxplot.{col}.pdf"), yscale=yscale, format_label=format_label,
                     prior=prior, var_name=var_name_mean)

    for (simu_model, rb_model), df_simu in trace_df.groupby(["simu", "model"]):
        plot_violin(df_simu, "Posterior", rename(f".violin.{simu_model}.{rb_model}.pdf"))
        plot_trace(df_simu, "Posterior", rename(f".trace.{simu_model}.{rb_model}.pdf"))

    output_df.to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
