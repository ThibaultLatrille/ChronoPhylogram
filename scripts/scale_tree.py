import argparse
from libraries import *


def parse_tree(tree_path: str) -> Tree:
    return rename_tree(prune_tree(open_tree(tree_path, format_ete3=1)))


def main(tree_1, tree_2, tree_output_1, tree_output_2):
    """
    Scale the tree_1 and tree_2 to have the same species and the same scale
    :return:
    """
    t_1 = parse_tree(tree_1)
    t_2 = parse_tree(tree_2)
    intersection = list(set(t_1.get_leaf_names()).intersection(set(t_2.get_leaf_names())))
    print(f"Intersection between tree_1 and tree_2 has {len(intersection)} leaves")
    pruned_1 = scale_tree(prune_tree(t_1, intersection))
    pruned_2 = scale_tree(prune_tree(t_2, intersection))
    assert len(pruned_1.get_leaf_names()) == len(
        pruned_2.get_leaf_names()), "Pruned trees do not have the same number of leaves"
    pruned_1.write(outfile=tree_output_1, format=1)
    pruned_2.write(outfile=tree_output_2, format=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--tree_1", help="Input tree 1", required=True)
    parser.add_argument("--tree_2", help="Input tree 2", required=True)
    parser.add_argument("--tree_output_1", help="Output tree 1", required=True)
    parser.add_argument("--tree_output_2", help="Output tree 2", required=True)
    args = parser.parse_args()
    main(args.tree_1, args.tree_2, args.tree_output_1, args.tree_output_2)
