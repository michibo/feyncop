
"""qcd_gen.py: This file is part of the feyncop/feyngen package.
    Implements functions to generate QCD graphs. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

import itertools

import phi_34_gen
from stuff import flip
from weighted_graph import WeightedGraph

fermion = 1
boson = 2


def gen_graphs(L, r_t2, u_t2, m, cntd, edge2cntd, vtx2cntd, notadpoles, chunk=None):
    """Generate QCD graphs with the desired parameters and properties.
        L: Loop number
        r_t2: Ext. fermion number
        u_t2: Ext. ghost number
        m: Ext. boson number"""

    phi34_graphs = phi_34_gen.gen_graphs(L, r_t2 + m + u_t2,
                                         cntd, edge2cntd, vtx2cntd,
                                         notadpoles, chunk=chunk)

    for g_phi34 in phi34_graphs:
        gen_qcd_graphs = (g.unlabeled_graph for g in gen_from_phi34_g(g_phi34, r_t2, u_t2, m))

        qcd_graphs = frozenset(gen_qcd_graphs)

        yield from qcd_graphs


def gen_from_phi34_g(fg, r_t2, u_t2, m):
    """Helper function: Generate full fledged QCD graphs from the bulk output of
        phi_34_gen.gen_graphs."""

    ext_vtcs = fg.external_vtcs_set
    int_vtcs = fg.internal_vtcs_set

    def dir_sign(e, v):
        v1, _ = fg.edges[e]
        return 1 if v1 == v else -1

    is_sl = [fg.is_selfloop(edge) for e, edge in enumerate(fg.edges)]
    ext_adj = [frozenset(fg.adj_edges(v, fg.edges_set)) for v in ext_vtcs]
    int_adj = [frozenset(fg.adj_edges(v, fg.edges_set)) for v in int_vtcs]
    for weights in itertools.product((1, 2), repeat=len(fg.edges)):
        fermion_edges = frozenset(e for e, w in enumerate(weights)
                                  if w == fermion)
        fermion_adj = [adj & fermion_edges for adj in int_adj]
        fermion_valences = (sum(2 if is_sl[e] else 1 for e in adj) for adj in fermion_adj)
        if any(val not in [0, 2] for val in fermion_valences):
            continue

        boson_edges = frozenset(e for e, w in enumerate(weights) if w == boson)
        boson_adj = [adj & boson_edges for adj in int_adj]
        boson_valences = (sum(2 if is_sl[e] else 1 for e in adj)
                          for adj in boson_adj)
        if any(val not in [1, 3, 4] for val in boson_valences):
            continue

        fermion_legs = sum(1 for adj in ext_adj for e in adj
                           if weights[e] == fermion)
        if fermion_legs != r_t2 + u_t2:
            continue

        boson_legs = sum(1 for adj in ext_adj for e in adj
                         if weights[e] == boson)
        if boson_legs != m:
            continue

        for fermion_weights in itertools.product((-1, 1), repeat=len(fermion_edges)):
            dir_weights = list(weights)
            for i, e in enumerate(fermion_edges):
                dir_weights[e] = fermion_weights[i]

            fermion_res = (sum(dir_sign(e, v) * dir_weights[e] for e in adj
                               if not is_sl[e])
                           for v, adj in zip(int_vtcs, fermion_adj))
            if any(fermion_res):
                continue

            edges = tuple(edge if w == 1 or w == 2 else flip(edge) for edge, w in zip(fg.edges, dir_weights))
            translated_weights = tuple(2 if w == 2 else 1 for w in weights)

            g = WeightedGraph(tuple(edges), tuple(translated_weights))

            fermion_loops = list(g.cntd_components_sub_edges(g.sub_edges_by_weight(fermion)))

            for ghost_loops in itertools.product((True, False), repeat=len(fermion_loops)):
                ghost_weights = list(g.edge_weights)
                for gh, loop in zip(ghost_loops, fermion_loops):
                    if gh:
                        for e in loop:
                            ghost_weights[e] = 3

                gw = WeightedGraph(tuple(edges), tuple(ghost_weights))
                if len(gw.external_vtcs_set & gw.vtcs_set_sub_edges(gw.sub_edges_by_weight(3) & gw.external_edges_set)) == u_t2:
                    yield gw
