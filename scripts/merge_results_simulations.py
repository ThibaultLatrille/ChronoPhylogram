import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(bayescode_list: str, output: str):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    list_df = []

    plt.figure(figsize=(12, 8))
    for path in bayescode_list:
        name_split = os.path.basename(path).split(".")[0].split("_")
        df = pd.read_csv(path.replace(".run", ".trace"), sep='\t')

        # Remove the first 50% of the rows
        plt.plot(df["lnprob"], label=name_split[1] + "_" + name_split[2])
        df = df.iloc[int(0.5 * len(df)):]
        for i, p in enumerate(["dataset", "gram", "seed"]):
            df[p] = name_split[i]
        list_df.append(df[["lnprob", "gram", "seed"]])
    df_out = pd.concat(list_df)
    plt.legend()
    plt.savefig(replace_last(output, ".pdf", "_trace.pdf"))
    df_out["dataset"] = df_out["gram"] + "_" + df_out["seed"]
    # Violin plot of the log-likelihood
    fig, axes = plt.subplots(2, 1, figsize=(16, 9))
    ax = axes[0]
    df_out.sort_values(by=["seed", "gram"], inplace=True)
    color_palette = {}
    for dataset, gram in zip(df_out["dataset"], df_out["gram"]):
        color_palette[dataset] = sns.color_palette("tab10")[0 if gram == "Chrono" else 1]
    sns.violinplot(x="dataset", y="lnprob", data=df_out, density_norm='width', palette=color_palette, ax=ax)
    # Add legend for the colors
    ax.plot([], [], color=sns.color_palette("tab10")[0], label="Chronogram")
    ax.plot([], [], color=sns.color_palette("tab10")[1], label="Phylogram")
    ax.legend()
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    ax = axes[1]
    for dataset, df_seed in df_out.groupby("seed"):
        phylo = df_seed[df_seed["gram"] == "Phylo"]["lnprob"].values
        chrono = df_seed[df_seed["gram"] == "Chrono"]["lnprob"].values
        assert len(phylo) == len(chrono), f"Length of phylo ({len(phylo)}) and chrono ({len(chrono)}) are different"
        diff = phylo - chrono
        sns.kdeplot(diff, label=dataset, ax=ax)
    ax.set_xlabel("lnprob(Phylo) - lnprob(Chrono)")
    # Vertical line at 0
    ax.axvline(0, color='black', linestyle='--')
    ax.legend(fontsize=6)
    plt.tight_layout()
    plt.savefig(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bayescode", nargs="+", help="Input BayesCode file", required=True)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.bayescode, args.output)
