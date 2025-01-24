import argparse
from gzip import open as gzopen
from ete3 import Tree


def open_tree(tree_path: str, format_ete3: int = 1) -> Tree:
    if tree_path.endswith(".gz"):
        newick = gzopen(tree_path).read().decode()
        return Tree(newick, format=format_ete3)
    else:
        return Tree(tree_path, format=format_ete3)


def prune_tree(input_tree: Tree, list_taxa: list = None) -> Tree:
    tree = input_tree.copy()

    # Prune tree
    if list_taxa is not None:
        tree.prune(list_taxa, preserve_branch_length=True)
        assert len(tree.get_leaves()) == len(list_taxa), f"Pruning failed: {len(tree.get_leaves())} != {len(list_taxa)}"

    # Add polytomies if branch length are 0
    remove_nodes = set([n for n in tree.traverse() if (n.dist == 0.0 and not n.is_root())])
    for n in remove_nodes:
        n.delete()
    assert len(set([n for n in tree.traverse()]).intersection(remove_nodes)) == 0, "Polytomies could not be removed"
    for n in tree.traverse():
        if not n.is_root():
            assert n.dist > 0.0, f"Branch length is 0.0 for node {n.name}"
    return tree


def main(input_path, output_path):
    # Open the tree
    input_tree = open_tree(input_path, format_ete3=1)
    # Prune the tree
    tree = prune_tree(input_tree)
    # Write the tree to a newick file
    tree.write(outfile=output_path, format=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input_tree", help="Input tree file", required=True)
    parser.add_argument("--output_tree", help="Output tree file", required=True)
    args = parser.parse_args()
    main(args.input_tree, args.output_tree)
