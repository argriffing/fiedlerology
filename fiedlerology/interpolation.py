"""
Interpolation on tree graphs.

This can be seen as a generalization of linear interpolation.

"""
from __future__ import division, print_function, absolute_import

import networkx as nx


def _get_first(it):
    for x in it:
        return x


def _weighted_mean(wv_pairs):
    wsum = 0
    wvsum = 0
    for w, v in wv_pairs:
        wsum += w
        wvsum += w*v
    return wvsum / wsum


def linalg_laplacian_schur_complement(G,
        observed_nodelist, unobserved_nodelist, weight='weight'):
    nodelist = observed_nodelist + unobserved_nodelist

def linalg_laplcian_interpolation(G,
        observed_nodelist, unobserved_nodelist,
        observed_values, weight='weight'):
    nodelist_keep = list(nodes)
    nodelist_remove = list(set(G) - set(nodes))
    nodelist = nodelist_keep + nodelist_remove
    L = nx.linalg.laplacian_matrix(G, nodelist).A
    k = len(nodelist_keep)
    v_kept = np.array(list(values), dtype=float)
    v_removed = -solve(L[k:, k:], L[k:, :k].dot(v_kept))
    return np.concatenate((v_kept, v_removed))


def arborescence_interpolation(A, root, v_in, weight='weight'):
    """
    Interpolation on an arborescence.

    A purely networkx based interface and implementation.

    Parameters
    ----------
    A : networkx DiGraph
        An arborescence.
    root : hashable
        The root of the arborescence.
    v_in : dict
        Sparse map from vertex to value.
        Keys not in G will be ignored.
    weight : string or None, optional
        The edge attribute that hold the numerical value used as a weight.

    Returns
    -------
    v_interp : dict
        Sparse map from vertex to value.
        All nodes in G should be represented as keys in this dict.

    """
    #TODO more input validation
    edges = list(nx.bfs_edges(A, root))

    # Do the pass from the leaves to the root.
    # For each non-root node, get the effective weight and value
    # of the subtree (including the edge leading to that node).
    eff_wv_pairs = {}
    for na, nb in reversed(edges):
        in_weight = A[na][nb][weight]
        if nb in v_in:
            eff_wv_pairs[nb] = (in_weight, v_in[nb])
        else:
            ncs = list(A[nb])
            eff_weight = None
            eff_value_total = None
            for nc in ncs:
                wv_pair = eff_wv_pairs[nc]
                if wv_pair is None:
                    continue
                w, v = wv_pair
                if eff_weight is None:
                    eff_weight = w
                    eff_value_total = w * v
                else:
                    eff_weight += w
                    eff_value_total += w * v
            if eff_weight is None:
                eff_wv_pairs[nb] = None
            else:
                eff_value = eff_value_total / eff_weight
                eff_wv_pairs[nb] = (eff_weight, eff_value)

    # Do the pass from the root to the leaves.
    v_interp = {}
    for na, nb in [(None, root)] + edges:
        if nb in v_in:
            v_interp[nb] = v_in[nb]
        else:
            wv_pairs = []
            if na is not None:
                w = A[na][nb][weight]
                v = v_interp[na]
                wv_pairs.append((w, v))
            for nc in A[nb]:
                wv_pair = eff_wv_pairs[nc]
                if wv_pair is not None:
                    wv_pairs.append(wv_pair)
            if not wv_pairs:
                raise Exception('failed to find an interpolation')
            v_interp[nb] = _weighted_mean(wv_pairs)

    return v_interp


def tree_interpolation(G, v, weight='weight'):
    """
    Interpolation on a tree domain.

    A purely networkx based interface and implementation.

    Parameters
    ----------
    G : networkx Graph
        Undirected weighted acyclic networkx graph.
    v : dict
        Sparse map from vertex to value.
        Keys not in G will be ignored.
    weight : string or None, optional
        The edge attribute that hold the numerical value used as a weight.

    Returns
    -------
    v_interp : dict
        Sparse map from vertex to value.
        All nodes in G should be represented as keys in this dict.

    """
    if weight is None:
        raise NotImplementedError(
                'interpolation is implemented only for weighted trees')
    nnodes = nx.number_of_nodes(G)
    if not nnodes:
        return {}
    if not nx.is_connected(G):
        raise ValueError('expected a tree, but the graph is not connected')
    edge_count = 0
    for na, nb in G.edges():
        if G[na][nb][weight] <= 0:
            raise ValueError('edge weights must be positive')
        edge_count += 1
    if edge_count != nnodes - 1:
        raise ValueError('expected a tree')
    observed_nodes = set(v) & set(G)
    # if no value is known on the tree then set all values to zero
    if not observed_nodes:
        return dict((node, 0) for node in G)
    vlocal = dict((node, v[node]) for node in observed_nodes)
    # If all values are known then no more work is needed.
    # Note that this will include the case of the single-node tree.
    if observed_nodes == set(G):
        return vlocal

    # construct the arborescence rooted at an observed node
    root = _get_first(observed_nodes)
    A = nx.bfs_tree(G, root)
    return arborescence_interpolation(A, root, vlocal, weight=weight)


def forest_interpolation(G, v, weight='weight'):
    """
    Interpolation on a forest domain.

    A purely networkx based interface and implementation.
    Separately interpolate each tree in the forest.

    Parameters
    ----------
    G : networkx Graph
        Undirected weighted networkx graph.
        This should be a forest (a graph whose connected components are trees).
    v : dict
        Sparse map from vertex to value.
    weight : string or None, optional
        The edge attribute that hold the numerical value used as a weight.

    Returns
    -------
    v_interp : dict
        Sparse map from vertex to value.

    """
    vlocal = {}
    for G_sub in nx.connected_component_subgraphs(G):
        vlocal.update(tree_interpolation(G, v, weight=weight))
    return vlocal

