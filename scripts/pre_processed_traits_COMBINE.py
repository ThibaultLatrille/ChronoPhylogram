import os
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
from libraries import *


def main(path_input_traits: str, path_input_tree: str, path_output_tree: str,
         path_output_traits: str, args_trait: str, sex: str):
    for path in [path_input_traits, path_input_tree]:
        assert os.path.exists(path), f"Path {path} does not exist"
    for path in [path_output_tree, path_output_traits]:
        os.makedirs(os.path.dirname(path), exist_ok=True)

    t = open_tree(path_input_tree, format_ete3=1)
    tree = name_internal_nodes(rename_tree(prune_tree(t)))
    set_taxa_names = set(tree.get_leaf_names())
    print(f"The tree has {len(set_taxa_names)} leaves")

    df_traits = pd.read_csv(path_input_traits)
    df_traits["Genus_Species"] = df_traits["iucn2020_binomial"].apply(lambda x: x.replace(" ", "_"))
    assert "Genus_Species" in df_traits.columns
    dico_trait_names = {"bodyMass": "adult_mass_g", "brainMass": "adult_brain_mass_g"}
    trait = dico_trait_names[args_trait]
    assert trait in df_traits.columns
    print(f"\nPhenotype considered is {trait}")
    print(f"The trait dataframe has {len(df_traits)} rows initially.")
    df_traits = df_traits[df_traits["Genus_Species"].isin(set_taxa_names)]
    print(f"The trait dataframe has {len(df_traits)} rows after filtering taxa also in the tree.")
    df_traits = df_traits[np.isfinite(df_traits[trait])]
    print(f"The trait dataframe has {len(df_traits)} rows after filtering not finite values.")
    # Log transform the trait
    df_traits[trait] = np.log(df_traits[trait])
    print(f"The trait is log transformed.")

    # Filter the tree and create dictionaries for variance and mean of traits
    set_taxa_names = set_taxa_names.intersection(set(df_traits["Genus_Species"]))
    print(f"\nThe intersection of the tree and trait dataframe has {len(set_taxa_names)} taxa.")
    tree = prune_tree(tree, list(set_taxa_names))
    assert len(tree.get_leaves()) == len(set_taxa_names)
    dico_traits = defaultdict(list)
    for taxa_name in set_taxa_names:
        dico_traits["TaxonName"].append(taxa_name)

    # Filter out the species with a unique row in the trait dataframe
    var_df = df_traits.copy()
    var_grouped = {k: v for k, v in var_df.groupby("Genus_Species")}
    print(f"The trait dataframe has {len(var_grouped)} taxa after keeping taxon with at least 1 individual.")
    for taxa_name in set_taxa_names:
        if taxa_name in var_grouped:
            phenotype_mean = np.average(var_grouped[taxa_name][trait])
            dico_traits[f"Phenotype_mean"].append(phenotype_mean)
    print(f"{len(var_grouped)} species with mean computed.")

    df_traits = pd.DataFrame(dico_traits)
    df_traits = df_traits[np.isfinite(df_traits.drop(["TaxonName"], axis=1)).any(axis=1)]
    df_traits.to_csv(path_output_traits, sep="\t", index=False, na_rep="NaN")

    # Prune the tree and write it
    set_taxa_names = set(df_traits["TaxonName"])
    tree = scale_tree(prune_tree(tree, list(set_taxa_names)))

    print(f"The final tree has {len(tree.get_leaves())} taxa.")
    tree_length = sum([node.dist for node in tree.traverse()])
    print(f"The tree length is {tree_length}.")
    tree.write(outfile=path_output_tree, format=3)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_traits", help="Input trait file", required=True)
    parser.add_argument("--input_tree", help="Input dS tree file", required=True)
    parser.add_argument("--output_tree", help="Output tree file", required=True)
    parser.add_argument("--output_traits", help="Output traits file", required=True)
    parser.add_argument('--trait', help="trait name", default="false", type=str, required=False)
    parser.add_argument('--sex', help="Sex (m or f)", default="m", type=str, required=False)
    args = parser.parse_args()
    main(args.input_traits, args.input_tree, args.output_tree, args.output_traits, args.trait, args.sex)
