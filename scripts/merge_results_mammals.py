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
        for i, p in enumerate(["dataset", "gram", "sex", "logT"]):
            df[p] = name_split[i]
        list_df.append(df[["lnprob", "gram", "sex"]])
    df_out = pd.concat(list_df)
    plt.legend()
    plt.savefig(replace_last(output, ".pdf", "_trace.pdf"))
    df_out["dataset"] = df_out["gram"] + "_" + df_out["sex"]
    # Violin plot of the log-likelihood
    plt.figure(figsize=(12, 8))
    sns.violinplot(x="dataset", y="lnprob", data=df_out)
    plt.savefig(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bayescode", nargs="+", help="Input BayesCode file", required=True)
    parser.add_argument("--output", help="Output pdf file", required=True)
    args = parser.parse_args()
    main(args.bayescode, args.output)
