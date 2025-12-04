import itertools
import numpy as np
from ete3 import Tree


def brownian_fitting(tree: Tree, trait: str) -> (float, float, float):
    leaves = tree.get_leaves()
    n = len(leaves)
    C = np.zeros((n, n))
    x = np.zeros(n)
    for i, leaf in enumerate(leaves):
        x[i] = float(getattr(leaf, trait))

    for node in tree.traverse("levelorder"):
        if node.is_root():
            node.t = 0
        else:
            node.t = node.dist + node.up.t

    for i, j in itertools.product(range(n), range(n)):
        if i == j:
            C[i, j] = leaves[i].t
        elif i < j:
            ancestor = leaves[i].get_common_ancestor(leaves[j])
            C[i, j] = ancestor.t
            C[j, i] = ancestor.t

    invC = np.linalg.inv(C)
    ones = np.ones(n)
    v = np.dot(ones.T, invC)
    anc_z = float(np.dot(v, x)) / float(np.dot(v, ones))
    assert np.isfinite(anc_z)

    d = (x - anc_z * ones)
    var = float(np.dot(d.T, np.dot(invC, d)) / (n - 1))

    ll = -0.5 * (n * np.log(2 * np.pi * var) + np.log(np.linalg.det(C)) + (1 / var) * np.dot(d.T, np.dot(invC, d)))

    return anc_z, var, float(ll)
