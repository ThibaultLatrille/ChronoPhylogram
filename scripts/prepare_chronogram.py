import os
import argparse
from libraries import *


def name_nodes(tree: Tree) -> Tree:
    node_i = 0
    for n in tree.traverse():
        if n.is_leaf():
            n.name = f"leaf_{n.name}"
        elif not n.is_leaf():
            n.name = f"node_{node_i}"
            node_i += 1
        if n.is_root():
            n.name = "Root"
            continue
        assert n.dist > 0.0
    # Check that the tree is well named
    assert len(set([n.name for n in tree.traverse()])) == len([n.name for n in tree.traverse()])
    return tree


def main(tree_input, tree_output):
    os.makedirs(os.path.dirname(tree_output), exist_ok=True)
    t = open_tree(tree_input, format_ete3=1)
    t = prune_tree(t)
    t = name_nodes(t)
    t.write(outfile=tree_output, format=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", help="Input tree", required=True)
    parser.add_argument("--output", help="Output tree", required=True)
    args = parser.parse_args()
    main(args.input, args.output)
