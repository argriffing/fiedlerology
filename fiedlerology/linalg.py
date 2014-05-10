"""
Interpolation on tree graphs.

This can be seen as a generalization of linear interpolation.

"""
from __future__ import division, print_function, absolute_import

import networkx as nx
import numpy as np
from scipy.linalg import solve, pinvh


def laplacian(G,
        observed_nodelist, unobserved_nodelist, weight='weight'):
    """
    """
    nodelist = observed_nodelist + unobserved_nodelist
    for nodes in observed_nodelist, unobserved_nodelist, nodelist:
        if len(set(nodes)) < len(nodes):
            raise ValueError('duplicate nodes in a list')
    if set(nodelist) != set(G):
        raise ValueError('node set mismatch')
    L = nx.linalg.laplacian_matrix(G, nodelist)
    return L.A


def laplacian_interpolator(G,
        observed_nodelist, unobserved_nodelist, weight='weight'):
    """
    """
    L = laplacian(G, observed_nodelist, unobserved_nodelist, weight=weight)
    k = len(observed_nodelist)
    return -solve(L[k:, k:], L[k:, :k])


def laplacian_interpolator_multiply(G,
        observed_nodelist, unobserved_nodelist,
        observed_values, weight='weight'):
    """
    """
    L = laplacian(G, observed_nodelist, unobserved_nodelist, weight=weight)
    k = len(observed_nodelist)
    v = np.array(list(values), dtype=float)
    return -solve(L[k:, k:], L[k:, :k].dot(v))


def laplacian_sc(G, observed_nodelist, unobserved_nodelist, weight='weight'):
    """
    Laplacian Schur complement.

    """
    L = laplacian(G, observed_nodelist, unobserved_nodelist, weight=weight)
    k = len(observed_nodelist)
    return L[:k, :k] - L[:k, k:].dot(solve(L[k:, k:], L[k:, :k]))


def laplacian_sc_multiply(G,
        observed_nodelist, unobserved_nodelist,
        observed_values, weight='weight'):
    """
    Action of the Schur complement.

    """
    L = laplacian(G, observed_nodelist, unobserved_nodelist, weight=weight)
    k = len(observed_nodelist)
    v = np.array(list(values), dtype=float)
    a = L[:k, :k].dot(v)
    b = L[:k, k:].dot(solve(L[k:, k:], L[k:, :k].dot(v)))
    return a - b


def laplacian_sc_pinv(G,
        observed_nodelist, unobserved_nodelist, weight='weight'):
    """
    Pseudo-inverse of Laplacian Schur complement.

    """
    sc = laplacian_schur_complement(G,
            observed_nodelist, unobserved_nodelist, weight=weight)
    return pinvh(sc)


def laplacian_sc_pinv_multiply(G,
        observed_nodelist, unobserved_nodelist,
        observed_values, weight='weight'):
    """
    Action of the pseudo-inverse of Laplacian Schur complement.

    """
    sc = laplacian_schur_complement(G,
            observed_nodelist, unobserved_nodelist, weight=weight)
    v = np.array(list(observed_values), dtype=float)
    x, residues, rank, s = lstsq(sc, v)
    return x

