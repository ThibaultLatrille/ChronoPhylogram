#!/usr/bin/env python3
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
    assert len(
        set([n for n in tree.traverse()]).intersection(remove_nodes)) == 0, "Some polytomies could not be removed"
    for n in tree.traverse():
        if not n.is_root():
            assert n.dist > 0.0, f"Branch length is 0.0 for node {n.name}"
    return tree


def rename_tree(tree: Tree) -> Tree:
    for node in tree.traverse("levelorder"):
        if node.is_leaf():
            node.name = node.name.replace(" ", "_")
    return tree


def name_internal_nodes(tree: Tree) -> Tree:
    node_i = 0
    for n in tree.traverse():
        if not n.is_leaf():
            n.name = f"node_{node_i}"
            node_i += 1
        if n.is_root():
            n.name = "Root"
            continue
        assert n.dist > 0.0
    return tree


def scale_tree(tree: Tree) -> Tree:
    d_total = sum([node.dist for node in tree.traverse()])
    for node in tree.traverse():
        node.dist = node.dist / d_total
    return tree
