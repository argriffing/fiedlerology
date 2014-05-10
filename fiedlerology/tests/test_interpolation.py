"""
"""
from __future__ import division, print_function, absolute_import

import networkx as nx
import numpy as np
from scipy.linalg import inv, solve, eigh
from scipy.sparse.linalg import eigsh
from numpy.testing import assert_allclose


import fiedlerology
from fiedlerology.interpolation import arborescence_interpolation


def _linalg_tree_interpolation(G, nodes, values):
    # only report values for the removed nodes
    nodelist_keep = list(nodes)
    nodelist_remove = list(set(G) - set(nodes))
    nodelist = nodelist_keep + nodelist_remove
    L = nx.linalg.laplacian_matrix(G, nodelist).A
    k = len(nodelist_keep)
    v_kept = np.array(list(values), dtype=float)
    v_removed = -solve(L[k:, k:], L[k:, :k].dot(v_kept))
    return np.concatenate((v_kept, v_removed))


def _linalg_tree_laplacian_multiply(G, nodes, values):
    nodelist_keep = list(nodes)
    nodelist_remove = list(set(G) - set(nodes))
    nodelist = nodelist_keep + nodelist_remove
    L = nx.linalg.laplacian_matrix(G, nodelist).A
    k = len(nodelist_keep)
    v_kept = np.array(list(values), dtype=float)
    L_schur = L[:k, :k] - L[:k, k:].dot(solve(L[k:, k:], L[k:, :k]))
    v_out = L_schur.dot(values)


def test_arborescence_interpolation():
    root = 2
    G = nx.Graph()
    G.add_weighted_edges_from([
        (root, 0, 1.0),
        (root, 1, 2.0),
        ])
    A = nx.bfs_tree(G, root)
    for na, nb in A.edges():
        A[na][nb]['weight'] = G[na][nb]['weight']
    v_in = {0 : 0.5, 1 : 0.6}
    v_out = arborescence_interpolation(A, root, v_in)
    print(v_out)

    # interpolation
    v_out_linalg = _linalg_tree_interpolation(G, (0, 1), (0.5, 0.6))
    print(v_out_linalg)

    #L = np.array([
        #[3.0, -1.0, -2.0],
        #[-1.0, 1.0, 0.0],
        #[-2.0, 0.0, 2.0],
        #], dtype=float)
    #v_foo = np.empty(3)
    #for node, value in v_out.items():
        #v_foo[node] = value
    #print(np.dot(L, v_foo))

    #v_out_linalg = _linalg_tree_interpolation(G, (1, 2), (0.5, 0.6))

