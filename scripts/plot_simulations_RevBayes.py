import os
import argparse
from glob import glob
import numpy as np
from natsort import natsorted
import pandas as pd
from collections import defaultdict
from os.path import basename, isdir
from libraries_plot import hist_plot


def open_log(path):
    df = pd.read_csv(path, sep="\t")
    if "alpha" in df.columns:
        df["t_half"] = np.log(2) / df["alpha"]
    if "var_rates" in df.columns and "sigma2_root" in df.columns:
        df["std_rates"] = np.sqrt(df["var_rates"] / df["sigma2_root"])
    return df


def posterior(df, col="alpha", burnin=0.5):
    return np.nanmean(df[col][int(len(df) * burnin):])


def main(folder, output):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    rb_models = ["simple_OU", "relaxed_BM"]
    model_prefs = {"moving_optimum": 0, "directional": 1, "neutral": 3}
    models_path = {basename(p): p for p in glob(folder + "/*") if isdir(p) and basename(p) in model_prefs}
    models = list(sorted(models_path, key=lambda x: model_prefs[x] if x in model_prefs else -1))

    replicates = {m: natsorted([i for i in glob(f"{p}/*") if isdir(i)]) for m, p in models_path.items()}
    assert len(set([len(g) for g in replicates.values()])) == 1

    parameters = {"num_rate_changes": ("N", "linear"), "std_rates": ("std rates", "log"),
                  "var_multiplier": ("var", "log"), "var_rates": ("var rates", "log"),
                  "is_BM": ("p[BM]", "uniform"), "is_OU": ("p[OU]", "uniform"),
                  "sigma2": ("sigma2", "log"), "theta": ("theta", "linear"),
                  "alpha": ("alpha", "log"), "t_half": ("t 1/2", "log")}

    trace_dict = defaultdict(lambda: defaultdict(list))
    for m in models:
        for f, folderpath in enumerate(replicates[m]):
            for rb in rb_models:
                trace_path = f"{folderpath}/{rb}_RJ.log.gz"
                if not os.path.exists(trace_path):
                    continue
                trace_df = open_log(trace_path)
                gram = os.path.basename(folderpath).split("_")[1]
                key_name = f"{m}_{gram}"
                trace_dict[f"model_{rb}"][key_name].append(m)
                trace_dict[f"gram_{rb}"][key_name].append(gram)
                for col in parameters.keys():
                    if col not in trace_df.columns:
                        continue
                    trace_dict[col][key_name].append(posterior(trace_df, col=col, burnin=0.5))

    rename = lambda x: output.replace(".tsv", x)
    for col, (x_label, xscale) in parameters.items():
        if col not in trace_dict:
            continue
        hist_plot(trace_dict[col], x_label, rename(f".{col}.pdf"), xscale=xscale)

    out_dict = defaultdict(list)
    for rb in rb_models:
        for m in trace_dict[f"model_{rb}"].keys():
            for col, values in trace_dict.items():
                out_dict[f"{col}"].extend(values[m])
            for col in parameters.keys():
                if col not in out_dict:
                    continue
                print(f"{m} {col}: {np.mean(out_dict[col])} (std:{np.std(out_dict[col])})")
    df = pd.DataFrame(out_dict)
    df.to_csv(output, sep="\t", index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder', required=True, type=str, dest="folder", help="Input folder")
    parser.add_argument('-o', '--output', required=True, type=str, dest="output", help="Output path")
    args = parser.parse_args()
    main(args.folder, args.output)
