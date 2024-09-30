
"""qed_gen.py: This file is part of the feyncop/feyngen package.
    Implements functions to generate QED graphs. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

from itertools import product

from weighted_graph import WeightedGraph
from stuff import flip
import phi_k_gen

fermion = 1
boson = 2


def gen_graphs(L, r_t2, m, cntd, edge2cntd, vtx2cntd, notadpoles, furry):
    """Generate QED graphs with the desired parameters and properties.
        L: Loop number
        r_t2: External fermion number
        m: External boson number"""

    phi3_graphs = (phi_k_gen.gen_graphs(L, 3, r_t2 + m,
                                        cntd, edge2cntd, vtx2cntd, notadpoles))

    for g_phi3 in phi3_graphs:
        gen_qed_graphs = (g.unlabeled_graph for g in gen_from_phi3_g(g_phi3, r_t2, m))

        qed_graphs = frozenset(gen_qed_graphs)

        for g in qed_graphs:
            if furry:
                _, cycles = g.cycle_decomposition(g.sub_edges_by_weight(fermion))
                if any(len(c) % 2 for c in cycles):
                    continue
            yield g


def gen_from_phi3_g(fg, ext_fermion, ext_boson):
    """Helper function: Generate full fledged QED graphs
    from the bulk output of phi_k_gen.gen_graphs."""

    allowed_valencies = {(2, 1)}  # QED triple vertex
    if ext_fermion:
        allowed_valencies.add((1, 0))
    if ext_boson:
        allowed_valencies.add((0, 1))

    ext_vtcs = fg.external_vtcs_set
    int_vtcs = fg.internal_vtcs_set

    def dir(e, v):
        v1, _ = fg.edges[e]
        return 1 if v1 == v else -1

    is_sl = [fg.is_selfloop(edge) for edge in fg.edges]
    edge_valence = [2 if is_sl[e] else 1 for e, edge in enumerate(fg.edges)]

    ext_adj = [frozenset(fg.adj_edges(v, fg.edges_set)) for v in ext_vtcs]
    int_adj = [frozenset(fg.adj_edges(v, fg.edges_set)) for v in int_vtcs]

    for weights in product((fermion, boson), repeat=len(fg.edges)):
        boson_edges = frozenset(e for e, w in enumerate(weights)
                                if w == boson)
        fermion_edges = frozenset(e for e, w in enumerate(weights)
                                  if w == fermion)
        adjacence = [(adj & fermion_edges, adj & boson_edges)
                     for adj in int_adj]

        # check the valences
        valences = ((sum(edge_valence[e] for e in adj_fermion),
                     sum(edge_valence[e] for e in adj_boson))
                    for adj_fermion, adj_boson in adjacence)
        if any(val not in allowed_valencies for val in valences):
            continue

        # check the legs
        if ext_fermion:
            fermion_legs = sum(1 for adj in ext_adj for e in adj
                               if weights[e] == fermion)
            if fermion_legs != ext_fermion:
                continue

        if ext_boson:
            boson_legs = sum(1 for adj in ext_adj for e in adj
                             if weights[e] == boson)
            if boson_legs != ext_boson:
                continue

        fermion_adj = [a for a, _ in adjacence]

        for fermion_weights in product((-1, 1), repeat=len(fermion_edges)):
            dir_weights = list(weights)
            for i, e in enumerate(fermion_edges):
                dir_weights[e] = fermion_weights[i]
            # weight -1 for reversed fermion arrow

            fermion_res = (sum(dir(e, v) * dir_weights[e] for e in adj
                               if not is_sl[e])
                           for v, adj in zip(int_vtcs, fermion_adj))
            if any(fermion_res):
                continue

            edges = tuple(edge if w == 1 or w == boson else flip(edge)
                          for edge, w in zip(fg.edges, dir_weights))
            translated_weights = tuple(2 if w == 2 else 1 for w in weights)

            yield WeightedGraph(tuple(edges), tuple(translated_weights))
