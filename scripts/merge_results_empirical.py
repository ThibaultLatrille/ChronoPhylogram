import os
import argparse
import pandas as pd


def replace_last(s: str, old: str, new: str) -> str:
    li = s.rsplit(old, 1)
    return new.join(li)


def main(tsv_list: str, output: str):
    os.makedirs(os.path.dirname(output), exist_ok=True)
    list_df = []
    sex_dico = {"m": "\\Male", "f": "\\Female", "b": "Mixed"}
    for path in tsv_list:
        df = pd.read_csv(path, sep='\t')
        df["dataset"] = os.path.basename(path).replace(".tsv", "").split("_")[-1]
        df["sex"] = df["simu"].apply(lambda x: sex_dico[x.split("_")[0]])
        df["trait"] = df["simu"].apply(lambda x: x.split("_")[1])
        list_df.append(df)
    df_out = pd.concat(list_df)
    # Sort by trait, dataset, sex, logT and then method
    df_out = df_out.sort_values(by=["dataset", "trait", "sex", "rb"], ascending=[False, True, True, False])
    # Fill missing values with NaN
    df_out = df_out.fillna("NaN")
    # Save the dataframe to a tsv file
    df_out.to_csv(output, sep="\t", index=False, float_format="%.3f")

    dico_label = {"dataset": "Dataset", "trait": "Trait", "sex": "Sex", "rb": "Method",
                  "ntax": "$n$", "is_nuc": "$\\pi$"}
    columns = [i for i in dico_label if i in df_out]
    df_out = df_out[columns]
    df_out = df_out.rename(columns={i: (dico_label[i] if i in dico_label else i) for i in columns})
    df_out.to_latex(replace_last(output, ".tsv", ".tex"), index=False, float_format="%.3g",
                    escape=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_tsv", nargs="+", help="Input tsv file", required=True)
    parser.add_argument("--output_tsv", help="Output tsv file", required=True)
    args = parser.parse_args()
    main(args.input_tsv, args.output_tsv)
