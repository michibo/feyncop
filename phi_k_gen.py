
"""phi_k_gen.py: This file is part of the feyncop/feyngen package.
    Implements functions to generate phi^k graphs. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

import itertools
from graph import Graph
import collections

import nauty_ctrl


def calc_gen_params(L, k, m, cntd, notadpoles):
    """Helper function: Calculate the parameters for the call of multig."""

    D = k - 2

    n = (m + 2 * (L-1))
    l = (k * (L-1) + m)

    if l % D != 0 or n % D != 0:
        return (0, 0, 0)

    n //= D
    l //= D
    bulk_n = n + m
    bulk_l = l + m

    min_edges = 0
    max_edges = bulk_l

    if notadpoles:
        min_edges = bulk_l
    elif cntd:
        min_edges = bulk_n - 1
    elif k % 2 == 1:
        min_edges = bulk_n // 2
    else:
        min_edges = (m + 1) // 2

    return bulk_n, min_edges, max_edges


def gen_graphs(L, k, m, cntd, edge2cntd, vtx2cntd, notadpoles):
    """Generate phi^k graphs with the desired parameters and properties.
        L: Loop number
        k: Valency
        m: Ext. legs"""

    cntd = cntd | edge2cntd | vtx2cntd
    edge2cntd = edge2cntd | vtx2cntd
    notadpoles = notadpoles | vtx2cntd

    bulk_n, min_edges, max_edges = calc_gen_params(L, k, m, cntd, notadpoles)

    if bulk_n <= 0:
        return

    for g_bulk in nauty_ctrl.gen_nauty_graphs(
        bulk_n, cntd, k, min_edges, max_edges):
        labeled_graphs = (g for g in gen_from_bulk_g(
                    g_bulk, frozenset(range(bulk_n)),
                    k, m, notadpoles))

        unlabeled_graphs = frozenset(g.unlabeled_graph for g in labeled_graphs)

        for g in unlabeled_graphs:
            if vtx2cntd and not g.is_vtx_2_connected:
                continue
            elif edge2cntd and not g.is_edge_2_connected:
                continue
            elif cntd and not g.is_connected:
                continue

            if not vtx2cntd and notadpoles and g.is_tadpole:
                continue

            yield g


def gen_from_bulk_g(g, vtcs_set, k, m, notadpoles):
    """Generate full fledged phi^k graphs from the bulk output of multig."""

    valences = g.valency_dict
    leg_candidates = frozenset(v for v in vtcs_set if valences[v] == 1)

    def gen_ext_vtcs():
        if len(leg_candidates) < m:
            return

        if len(leg_candidates) == m:
            yield leg_candidates
            return

        if not notadpoles and k % 2 == 1:  # Extra legs can still be closed with self-loops!
            for ext_vtcs in itertools.combinations(leg_candidates, m):
                yield frozenset(ext_vtcs)

    for ext_vtcs in gen_ext_vtcs():
        int_vtcs = vtcs_set - ext_vtcs

        degree_defs = [k - valences[v] for v in int_vtcs]
        if any(d % 2 != 0 or d < 0 for d in degree_defs):
            continue

        selfloop_edges = [(v,v) for v,d in zip(int_vtcs, degree_defs) for i in range(d//2)]

        edges = g.edges + selfloop_edges
        yield Graph(edges)
