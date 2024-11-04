import argparse
import numpy as np
from collections import defaultdict
import pandas as pd
from libraries import open_tree


def clean_newick(newick: str) -> str:
    newick = newick.replace(",traits=", ":traits=")
    while "[&index" in newick:
        # find the first occurrence of "[&&NHX"
        start = newick.find("[&index")
        # find the corresponding closing bracket
        end = start + 7
        while newick[end] != "]":
            end += 1
        # remove the "[&&NHX...]" block
        assert newick[end + 1] == ":", "Closing bracket not found"
        pre = newick[:start]
        NHX = newick[start:end + 1].replace("[&index", "[&&NHX:index").replace(",", ":")
        assert NHX.startswith("[&&NHX")
        assert NHX.endswith("]")
        post = newick[end + 1:]
        # Find the branch length
        start = post.find(":")
        end = start + 1
        while post[end] not in [",", ")", ":", "[", "]", ";"]:
            end += 1
        branch_length = post[start:end]
        assert branch_length.startswith(":"), f"Branch length does not start with ':': {branch_length}"
        newick = pre + branch_length + NHX + post[end:]
    return newick


def main(input_trace: str, output_path: str):
    df = pd.read_csv(input_trace, sep="\t")
    dict_d = defaultdict(list)
    print(f"The trace contains {len(df)} samples")
    for tree in df["Tree"]:
        d_tree = open_tree(clean_newick(tree), format_ete3=1)
        for i, node in enumerate(d_tree.traverse()):
            node_n = node.name if node.is_leaf() else "NODE_" + str(i)
            dict_d[node_n].append(float(node.traits))
    # Copy only branch lengths and node names
    tree = open_tree(clean_newick(df["Tree"][0]), format_ete3=1)
    for i, node in enumerate(tree.traverse()):
        # Mean, 5% and 95% quantiles
        node_n = node.name if node.is_leaf() else "NODE_" + str(i)
        node.add_feature("Phenotype_mean", np.mean(dict_d[node_n], axis=0))
        node.add_feature("Phenotype_mean_min", np.quantile(dict_d[node_n], 0.05, axis=0))
        node.add_feature("Phenotype_mean_max", np.quantile(dict_d[node_n], 0.95, axis=0))

    tree.write(outfile=output_path, format=1, features=["Phenotype_mean", "Phenotype_mean_min", "Phenotype_mean_max"])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input trace file", required=True)
    parser.add_argument("--output", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.input.replace(".log.gz", ".trees.gz"), args.output)
