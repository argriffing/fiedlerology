"""
Interpolation on tree graphs.

This can be seen as a generalization of linear interpolation.

"""
from __future__ import division, print_function, absolute_import


def get_first(it):
    for x in it:
        return x


def nx_arborescence_interpolate(A, root, v, weight=None):
    """
    Interpolation on an arborescence.

    A purely networkx based interface and implementation.

    Parameters
    ----------
    A : networkx DiGraph
        An arborescence.
    root : hashable
        The root of the arborescence.
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
    # Do the pass from the leaves to the root.
    # For each edge, get the effective weight and value.
    edge_eff_weight = {}
    edge_eff_value = {}
    for 

    # do the pass from the root to the leaves
    pass


def nx_tree_interpolate(G, v, weight=None):
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
    root = get_first(observed_nodes)
    A = nx.bfs_tree(G, root)
    return nx_arborescence_interpolate(A, root, vlocal, weight=weight)


def nx_forest_interpolate(G, v, weight=None):
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
        vlocal.update(nx_tree_interpolate(G, v, weight=weight))
    return vlocal

