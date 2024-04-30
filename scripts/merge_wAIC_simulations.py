import os
import argparse
from collections import defaultdict
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def parse_bayescode_list(bayescode_list: list[str], branch=False) -> pd.DataFrame:
    list_df = []
    for path in bayescode_list:
        name_split = os.path.basename(path).split(".")[0].split("_")
        df = pd.read_csv(path, sep='\t')
        df = df[(df["name"] != "wAIC") if branch else (df["name"] == "wAIC")]
        for i, p in enumerate(["dataset", "gram", "seed"]):
            df[p] = name_split[i]
        list_df.append(df[["name", "Exp", "Var", "gram", "seed"]])
    df_out = pd.concat(list_df)
    df_out["dataset"] = df_out["gram"] + "_" + df_out["seed"]
    df_out.sort_values(by=["seed", "gram"], inplace=True)
    df_dict = defaultdict(list)
    for seed, df_seed in df_out.groupby("seed"):
        phylo = df_seed[df_seed["gram"] == "Phylo"]
        chrono = df_seed[df_seed["gram"] == "Chrono"]
        assert len(phylo) == len(chrono), f"Phylo and Chrono have different number of rows for seed {seed}"
        delta_waic = chrono["Exp"] - phylo["Exp"]
        if branch:
            delta_waic -= chrono["Var"] - phylo["Var"]
        df_dict["ΔwAIC"].extend(delta_waic)
        df_dict["experiment"].extend([seed] * len(delta_waic))
        df_dict["name"].extend(chrono["name"])
    # Plot delta wAIC for each experiment
    return pd.DataFrame(df_dict)


def main(bayescode_list: list[str], output: str):
    os.makedirs(os.path.dirname(output), exist_ok=True)

    df = parse_bayescode_list(bayescode_list, branch=False)
    # Violin plot of the log-likelihood
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    sns.violinplot(data=df, x="experiment", y="ΔwAIC", inner="quart", fill=True, ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_xlabel("Experiment")
    ax.set_ylabel("ΔwAIC (Chrono - Phylo)")
    ax.axhline(0, color="black", linestyle="--")
    ax.axhline(1, color="red", linestyle="--")
    ax.axhline(-1, color="red", linestyle="--")
    plt.tight_layout()
    ax.set_title("wAIC comparison")
    plt.tight_layout()
    plt.savefig(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bayescode", nargs="+", help="Input BayesCode file", required=True)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.bayescode, args.output)
