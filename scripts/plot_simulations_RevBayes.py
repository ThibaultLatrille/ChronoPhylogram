import os
import argparse
from glob import glob
import numpy as np
from natsort import natsorted
import pandas as pd
from collections import defaultdict
from os.path import basename, isdir
from libraries_plot import vert_boxplot
from merge_results_simulations import plot_violin, plt, my_dpi


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
                for col in trace_df.columns:
                    if col not in parameters:
                        continue
                    post_dict[col][key_name].append(posterior(trace_df, col=col))

    df_out = pd.concat(list_df)
    return df_out, post_dict


def plot_trace(df_out: pd.DataFrame, col: str, output: str):
    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    for _, df in df_out.groupby(["seed", "gram"]):
        plt.plot(df[col])
    plt.savefig(output)
    plt.clf()


def format_label(l):
    rm_list = {"Both", "SBM", "RBM", "Switch", "RJ", "SOU", "ROU", "neutral", "moving", "multi", "optimum"}
    for rm in rm_list:
        l = l.replace(rm, "")
    while "  " in l:
        l = l.replace("  ", " ")
    return l.replace("Chrono", "Chronogram").replace("Phylo", "Phylogram")


def main(folder, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    rb_models = ["relaxed_BM_RJ", "relaxed_OU_RJ", "simple_OU_RJ", "simple_BM_REML", "simple_BM_nodes",
                 "simple_BM_SwitchREML", "simple_BM_Switchnodes"]
    simu_model_prefs = {"moving_optimum": 0, "multi_optimum": 1, "neutral": 3}
    simu_models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p) and basename(p) in simu_model_prefs}
    simu_models = list(sorted(simu_models_path, key=lambda x: simu_model_prefs[x] if x in simu_model_prefs else -1))

    replicates = {m: natsorted([i for i in glob(f"{p}/*") if isdir(i)]) for m, p in simu_models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    parameters = {"num_rate_changes": ("N", "linear"), "std_rates": ("std rates", "log"),
                  "var_multiplier": ("var", "log"), "var_rates": ("var rates", "log"),
                  "num_theta_changes": ("N", "linear"), "theta_multiplier": ("theta", "log"),
                  "is_BM": ("Probability of Brownian\n(not Ornstein-Uhlenbeck)", "uniform"),
                  "is_OU": ("Probability of OU\n(not Brownian)", "uniform"),
                  "is_nuc": ("Probability of Phylogram", "uniform"),
                  "sigma": ("sigma", "log"), "theta": ("theta", "linear"),
                  "alpha": ("alpha", "log"), "t_half": ("t 1/2", "log")}
    trace_df, post_dict = get_dicts(simu_models, replicates, rb_models, parameters)

    rename = lambda x: output.replace(".tsv", x)
    for col, dict_input in post_dict.items():
        y_label, yscale = parameters[col]
        vert_boxplot(dict_input, y_label, rename(f".boxplot.{col}.pdf"), yscale=yscale, format_label=format_label)

    for (simu_model, rb_model), df_simu in trace_df.groupby(["simu", "model"]):
        plot_violin(df_simu, "Posterior", rename(f".violin.{simu_model}.{rb_model}.pdf"))
        plot_trace(df_simu, "Posterior", rename(f".trace.{simu_model}.{rb_model}.pdf"))
    trace_df.to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
