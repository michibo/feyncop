
"""yukawa_phi4_gen.py: This file is part of the feyncop/feyngen package.
    Implements functions to generate Yukawa+Phi4 graphs. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

from itertools import product
from weighted_graph import WeightedGraph

import phi_34_gen

# 1 stands for fermion, 2 for boson
fermion = 1
boson = 2


def gen_graphs_yukawa_phi4(loops, ext_fermion, ext_boson,
                           cntd, edge2cntd, vtx2cntd, notadpoles):
    """
    Generate Yukawa+Phi^4 graphs with the desired parameters and properties.

    ``loops``: loop number
    ``ext_fermion``: External fermion number
    ``ext_boson``: External boson number

    EXAMPLES::

        sage: from yukawa_phi4_gen import *
        sage: L = gen_graphs_yukawa_phi4(3,0,0,True,False,False,True)
        sage: list(L)
        [G[[1,0,A],[1,0,A],[1,0,A],[1,0,A]]/48]
    """
    phi34_graphs = phi_34_gen.gen_graphs(loops, ext_fermion + ext_boson,
                                         cntd, edge2cntd, vtx2cntd, notadpoles)

    for g_phi34 in phi34_graphs:
        gen_qcd_graphs = (g.unlabeled_graph
                          for g in gen_yukawa_phi4_from_phi34(g_phi34,
                                                              ext_fermion, ext_boson))
        yield from frozenset(gen_qcd_graphs)


def gen_yukawa_phi4_from_phi34(graph, ext_fermion, ext_boson):
    """
    Helper function: Generate full fledged Yukawa-Phi4 graphs from the bulk output of
    phi_34_gen.gen_graphs.

    graph: phi34 graph
    ext_fermion: external fermion number
    ext_boson: external boson number
    """
    ext_vtcs = graph.external_vtcs_set
    int_vtcs = graph.internal_vtcs_set

    is_sl = [graph.is_selfloop(edge) for edge in graph.edges]
    edge_valence = [2 if is_sl[e] else 1 for e, edge in enumerate(graph.edges)]

    ext_adj = [frozenset(graph.adj_edges(v, graph.edges_set)) for v in ext_vtcs]
    int_adj = [frozenset(graph.adj_edges(v, graph.edges_set)) for v in int_vtcs]

    allowed_valencies = {(0, 4), (2, 1)}  # phi^4 and Yukawa
    if ext_fermion:
        allowed_valencies.add((1, 0))
    if ext_boson:
        allowed_valencies.add((0, 1))

    for weights in product((fermion, boson), repeat=len(graph.edges)):
        fermion_edges = frozenset(e for e, w in enumerate(weights)
                                  if w == fermion)
        boson_edges = frozenset(e for e, w in enumerate(weights)
                                if w == boson)
        adjacence = [(adj & fermion_edges, adj & boson_edges)
                     for adj in int_adj]
        valences = ((sum(edge_valence[e] for e in adj_fermion),
                     sum(edge_valence[e] for e in adj_boson))
                    for adj_fermion, adj_boson in adjacence)
        if any(val not in allowed_valencies for val in valences):
            continue

        fermion_legs = sum(1 for adj in ext_adj for e in adj if weights[e] == fermion)
        if fermion_legs != ext_fermion:
            continue

        boson_legs = sum(1 for adj in ext_adj for e in adj if weights[e] == boson)
        if boson_legs != ext_boson:
            continue

        fermion_adj = [a for a, _ in adjacence]
        for fermion_weights in product((-1, 1), repeat=len(fermion_edges)):
            # handling the directions
            dir_weights = list(weights)
            for i, e in enumerate(fermion_edges):
                dir_weights[e] = fermion_weights[i]

            def dir(e, v):
                v1, v2 = graph.edges[e]
                return 1 if v1 == v else -1

            fermion_res = (sum(dir(e, v) * dir_weights[e] for e in adj
                               if not is_sl[e])
                           for v, adj in zip(int_vtcs, fermion_adj))
            if any(fermion_res):
                continue

            def flip(xy):
                x, y = xy
                return (y, x)

            edges = tuple(edge if w == 1 or w == 2 else flip(edge)
                          for edge, w in zip(graph.edges, dir_weights))
            translated_weights = tuple(2 if w == 2 else 1 for w in weights)

            g = WeightedGraph(edges, translated_weights)

            _, cycles = g.cycle_decomposition(g.sub_edges_by_weight(fermion))
            if any(len(c) % 2 for c in cycles):
                continue
            yield g
