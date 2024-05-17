import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

my_dpi = 150
fontsize = 24
fontsize_legend = 18


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(bayescode_list: str, output: str):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    list_df = []

    plt.figure(figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    col = "multiLnprob"
    for path in bayescode_list:
        name_split = os.path.basename(path).split(".")[0].split("_")
        df = pd.read_csv(path.replace(".run", ".trace.gz"), sep='\t')

        # Remove the first 50% of the rows
        plt.plot(df[col], label=name_split[1] + "_" + name_split[2])
        df = df.iloc[int(0.5 * len(df)):]
        for i, p in enumerate(["dataset", "gram", "seed"]):
            df[p] = name_split[i]
        list_df.append(df[[col, "gram", "seed"]])
    df_out = pd.concat(list_df)
    plt.legend(fontsize=fontsize_legend)
    plt.savefig(replace_last(output, ".pdf", "_trace.pdf"))
    df_out["dataset"] = df_out["gram"] + "_" + df_out["seed"]
    # Violin plot of the log-likelihood
    fig, axes = plt.subplots(2, 1, figsize=(1920 / my_dpi, 1080 / my_dpi), dpi=my_dpi)
    ax = axes[0]
    df_out.sort_values(by=["seed", "gram"], inplace=True)
    palette = {"Chrono": sns.color_palette("tab10")[0], "Phylo": sns.color_palette("tab10")[1]}
    sns.violinplot(data=df_out, x="seed", y=col, hue="gram", split=True, inner="quart", fill=True,
                   density_norm='width', palette=palette, ax=ax)
    ax.legend(fontsize=fontsize_legend)
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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bayescode", nargs="+", help="Input BayesCode file", required=True)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.bayescode, args.output)
