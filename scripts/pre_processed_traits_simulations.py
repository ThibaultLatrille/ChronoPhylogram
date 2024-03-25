import argparse
from collections import defaultdict
import pandas as pd
from libraries import open_tree


def main(input_path, traits_path):
    tree = open_tree(input_path, format_ete3=1)

    dico_traits = defaultdict(list)
    for n in tree.get_leaves():
        dico_traits["TaxonName"].append(n.name)
        pheno_mean = float(getattr(n, "Phenotype_mean"))
        dico_traits["Phenotype_mean"].append(pheno_mean if pheno_mean != 0.0 else "NaN")
    df_traits = pd.DataFrame(dico_traits)
    df_traits.to_csv(traits_path, sep=("\t" if traits_path.endswith(".tsv") else ","), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tree file", required=True)
    parser.add_argument("--traitsfile", help="Output traits file", required=True)
    args = parser.parse_args()
    main(args.input, args.traitsfile)
