
"""phi_34_gen.py: This file is part of the feyncop/feyngen package.
    Implements functions to generate phi^3+phi^4 graphs. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

import itertools
from graph import Graph

import nauty_ctrl


def calc_gen_params(L, m, cntd, notadpoles):
    """Helper function: Calculate the parameters for the call of multig."""

    min_n = (m + 2 * (L - 1) + 1) // 2 + m
    max_n = (m + 2 * (L - 1)) + m

    min_l = (4 * (L - 1) + m + 1) // 2 + m // 2
    max_l = 3 * (L - 1) + 2 * m

    if not notadpoles:
        if cntd:
            min_l = min_n - 1
        else:
            min_l = (m + 1) // 2

    return min_n, max_n, min_l, max_l


def gen_graphs(L, m, cntd, edge2cntd, vtx2cntd, notadpoles):
    """Generate phi^3 + phi^4 graphs with the desired parameters and
        properties.
        L: Loop number
        k: Valency
        m: Ext. legs"""

    cntd = cntd | edge2cntd | vtx2cntd
    edge2cntd = edge2cntd | vtx2cntd
    notadpoles = notadpoles | vtx2cntd

    min_n, max_n, min_edges, max_edges = calc_gen_params(L, m, cntd, notadpoles)

    if max_n <= 0:
        return

    for bulk_n in range(min_n, max_n + 1):
        for g_bulk in nauty_ctrl.gen_nauty_graphs(
                bulk_n, cntd, 4, min_edges, max_edges):
            labeled_graphs = (g for g in gen_from_bulk_g(
                g_bulk, frozenset(range(bulk_n)),
                L, m, notadpoles))

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


def gen_from_bulk_g(g, vtcs_set, L, m, notadpoles):
    """Generate full fledged phi^3 + phi^4 graphs from the bulk output of
        multig."""

    valences = g.valency_dict
    leg_candidates = frozenset(v for v in vtcs_set if valences[v] == 1)

    def gen_ext_vtcs():
        if len(leg_candidates) < m:
            return

        if len(leg_candidates) == m:
            yield leg_candidates
            return

        if not notadpoles:  # Extra legs can still be closed with self-loops!
            for ext_vtcs in itertools.combinations(leg_candidates, m):
                yield frozenset(ext_vtcs)

    for ext_vtcs in gen_ext_vtcs():
        int_vtcs = vtcs_set - ext_vtcs

        degree_defs = [4 - valences[v] for v in int_vtcs]
        if any(d < 0 for d in degree_defs):
            continue

        selfloop_edges = [(v, v) for v, d in zip(int_vtcs, degree_defs)
                          for i in range(d // 2)]

        edges = g.edges + selfloop_edges
        if len(edges) - len(vtcs_set) == L - 1:
            yield Graph(edges)
