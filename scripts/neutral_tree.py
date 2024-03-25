import argparse
from libraries import *


def main(input_neutral_path, tree_path):
    neutral_tree = open_tree(input_neutral_path, format_ete3=1)
    for n in neutral_tree.traverse():
        n.dist = 0. if n.is_root() else float(getattr(n, "d"))
    neutral_tree = prune_tree(neutral_tree)
    # Write the topology to a newick file
    neutral_tree.write(outfile=tree_path, format=3)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--neutral_tree", help="Input neutral tree file", required=True)
    parser.add_argument("--tree", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.neutral_tree, args.tree)
