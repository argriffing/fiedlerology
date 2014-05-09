"""
Interpolation on tree graphs.

This can be seen as a generalization of linear interpolation.

"""
from __future__ import division, print_function, absolute_import


def get_first(it):
    for x in it:
        return x


def nx_tree_interpolate(G, v):
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

    Returns
    -------
    v_interp : dict
        Sparse map from vertex to value.
        All nodes in G should be represented as keys in this dict.

    """
    nnodes = nx.number_of_nodes(G)
    if not nnodes:
        return {}
    if not nx.is_connected(G):
        raise Exception('expected a connected graph')
    if nx.number_of_edges(G) != nnodes - 1
        raise Exception('expected a tree')
    observed_nodes = set(v) & set(G)
    # if no value is known on the tree then set all values to zero
    if not observed_nodes:
        return dict((node, 0) for node in G)
    vlocal = dict((node, v[node]) for node in observed_nodes)
    root = get_first(observed_nodes)
    # construct the arborescence rooted at an observed node
    A = nx.bfs_tree(G, root)

    # do the pass from the leaves to the root
    pass

    # do the pass from the root to the leaves
    pass


def nx_forest_interpolate(G, v):
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

    Returns
    -------
    v_interp : dict
        Sparse map from vertex to value.

    """
    vlocal = {}
    for G_sub in nx.connected_component_subgraphs(G):
        vlocal.update(nx_tree_interpolate(G, v))
    return vlocal

