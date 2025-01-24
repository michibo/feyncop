
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
        m: Ext. legs"""

    cntd = cntd | edge2cntd | vtx2cntd
    edge2cntd = edge2cntd | vtx2cntd
    notadpoles = notadpoles | vtx2cntd

    min_n, max_n, min_edges, max_edges = calc_gen_params(L, m,
                                                         cntd, notadpoles)

    if max_n <= 0:
        return

    for bulk_n in range(min_n, max_n + 1):
        for g_bulk in nauty_ctrl.gen_nauty_graphs(
                bulk_n, cntd, 4, min_edges, max_edges):
            labeled_graphs = (g for g in gen_from_bulk_g(
                g_bulk, frozenset(range(bulk_n)),
                L, m, notadpoles))

            # store graphs according to number of edges
            stock = {k: set() for k in range(min_edges, max_edges + 1)}

            for nedges, g in labeled_graphs:
                if vtx2cntd and not g.is_vtx_2_connected:
                    continue
                if edge2cntd and not g.is_edge_2_connected:
                    continue
                # useless, as we already asked nauty for connected graphs
                # elif cntd and not g.is_connected:
                #    continue
                if not vtx2cntd and notadpoles and g.is_tadpole:
                    continue
                new_g = g.unlabeled_graph
                if new_g not in stock[nedges]:
                    stock[nedges].add(new_g)
                    yield new_g


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

        if not notadpoles:
            # Extra legs can still be closed by adding self-loops!
            for ext_vtcs in itertools.combinations(leg_candidates, m):
                yield frozenset(ext_vtcs)

    for ext_vtcs in gen_ext_vtcs():
        int_vtcs = vtcs_set - ext_vtcs

        degree_defs = [4 - valences[v] for v in int_vtcs]
        if any(d < 0 for d in degree_defs):
            continue

        # adding loops to inner vertices of valence 1 or 2 (default 3 or 2)
        selfloop_edges = [(v, v) for v, d in zip(int_vtcs, degree_defs)
                          for i in range(d // 2)]

        edges = g.edges + selfloop_edges
        if len(edges) - len(vtcs_set) == L - 1:
            yield (len(edges), Graph(edges))


# sage versions


def gen_graphs_sage(L, m, cntd, edge2cntd, vtx2cntd, notadpoles):
    """Generate phi^3 + phi^4 graphs with the desired parameters and
        properties.
        L: Loop number
        m: Ext. legs"""

    cntd = cntd | edge2cntd | vtx2cntd
    edge2cntd = edge2cntd | vtx2cntd
    notadpoles = notadpoles | vtx2cntd

    min_n, max_n, min_edges, max_edges = calc_gen_params(L, m,
                                                         cntd, notadpoles)

    if max_n <= 0:
        return

    for bulk_n in range(min_n, max_n + 1):
        for g_bulk in nauty_ctrl.gen_nauty_graphs_sage(
                bulk_n, cntd, 4, min_edges, max_edges):
            labeled_graphs = (g for g in gen_from_bulk_g_sage(
                g_bulk, frozenset(range(bulk_n)),
                L, m, notadpoles))

            # store graphs according to number of edges
            stock = {k: set() for k in range(min_edges, max_edges + 1)}

            for nedges, g in labeled_graphs:
                if vtx2cntd and not g.is_biconnected():
                    continue

                # TODO
                # if edge2cntd and g.edge_connectivity(value_only=True) < 2:
                #    continue

                # useless, as we already asked nauty for connected graphs
                # elif cntd and not g.is_connected:
                #    continue

                # TODO:
                # if not vtx2cntd and notadpoles and g.is_tadpole:
                #    continue

                new_g = g.canonical_label().copy(immutable=True)
                if new_g not in stock[nedges]:
                    stock[nedges].add(new_g)
                    yield new_g


def gen_from_bulk_g_sage(g, vtcs_set, L, m, notadpoles):
    leg_candidates = frozenset(v for v in g if g.degree(v) == 1)

    def gen_ext_vtcs():
        if len(leg_candidates) < m:
            return

        if len(leg_candidates) == m:
            yield leg_candidates
            return

        if not notadpoles:
            # Extra legs can still be closed by adding self-loops!
            for ext_vtcs in itertools.combinations(leg_candidates, m):
                yield frozenset(ext_vtcs)

    for ext_vtcs in gen_ext_vtcs():

        degree_defs = [(v, 4 - g.degree(v)) for v in g if v not in ext_vtcs]
        if any(d < 0 for d, _ in degree_defs):
            continue

        # adding loops to inner vertices of valence 1 or 2 (default 3 or 2)
        selfloop_edges = [(v, v) for v, d in degree_defs
                          for i in range(d // 2)]

        num_edges = g.num_edges() + len(selfloop_edges)
        if num_edges - g.num_verts() == L - 1:
            new_g = g.copy()
            new_g.add_edges(selfloop_edges)
            yield (num_edges, new_g)
