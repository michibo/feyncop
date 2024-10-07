
"""yukawa_phi4_gen.py: This file is part of the feyncop/feyngen package.
    Implements functions to generate Yukawa+Phi4 graphs. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

from itertools import product

from weighted_graph import WeightedGraph
from stuff import flip
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

        sage: L = gen_graphs_yukawa_phi4(0,0,4,True,False,False,True)
        sage: list(L)
        [G[[1,0,A],[2,0,A],[3,0,A],[4,0,A]]/24]
    """
    phi34_graphs = phi_34_gen.gen_graphs(loops, ext_fermion + ext_boson,
                                         cntd, edge2cntd, vtx2cntd, notadpoles)

    for g_phi34 in phi34_graphs:
        gen_yukawa_graphs = (g.unlabeled_graph
                             for g in gen_yukawa_phi4_from_phi34(g_phi34,
                                                                 ext_fermion,
                                                                 ext_boson))
        yield from frozenset(gen_yukawa_graphs)


def gen_yukawa_phi4_from_phi34(graph, ext_fermion, ext_boson):
    """
    Helper function: Generate full fledged Yukawa-Phi4 graphs
    from the bulk output of phi_34_gen.gen_graphs.

    ``graph``: phi34 graph
    ``ext_fermion``: external fermion number
    ``ext_boson``: external boson number

    EXAMPLES::

        sage: from yukawa_phi4_gen import *
        sage: from graph import Graph
        sage: G = Graph([[0,1],[0,1],[0,1],[0,1]])

        sage: L = gen_yukawa_phi4_from_phi34(G, 0, 0)
        sage: list(L)
        [G[[0,1,A],[0,1,A],[0,1,A],[0,1,A]]]

    """
    allowed_valencies = {(0, 4), (2, 1)}  # phi^4 and Yukawa
    if ext_fermion:
        allowed_valencies.add((1, 0))
    if ext_boson:
        allowed_valencies.add((0, 1))

    ext_vtcs = graph.external_vtcs_set
    int_vtcs = graph.internal_vtcs_set
    phi4_vtcs = [v for v in int_vtcs
                 if graph.vtx_valence(v, graph.edges_set) == 4]

    is_sl = [graph.is_selfloop(edge) for edge in graph.edges]
    edge_valence = [2 if is_sl[e] else 1 for e, edge in enumerate(graph.edges)]

    ext_adj = [frozenset(graph.adj_edges(v, graph.edges_set))
               for v in ext_vtcs]
    int_adj = [frozenset(graph.adj_edges(v, graph.edges_set))
               for v in int_vtcs]

    forced_edges = frozenset(i for i, e in enumerate(graph.edges)
                             if e[0] in phi4_vtcs or e[1] in phi4_vtcs)
    other_edges = [e for e in graph.edges_set if e not in forced_edges]

    def dir_sign(e, v):
        v1, _ = graph.edges[e]
        return 1 if v1 == v else -1

    for weights in product((fermion, boson), repeat=len(other_edges)):
        fermion_edges = frozenset(e for e, w in zip(other_edges, weights)
                                  if w == fermion)
        boson_edges = {e for e, w in zip(other_edges, weights)
                       if w == boson}
        boson_edges.update(forced_edges)

        full_weights = [boson] * len(graph.edges_set)
        pos = 0
        for i in graph.edges_set:
            if i in other_edges:
                full_weights[i] = weights[pos]
                pos += 1

        # check that vertices are correct
        adjacence = [(adj & fermion_edges, adj & boson_edges)
                     for adj in int_adj]
        valences = ((sum(edge_valence[e] for e in adj_fermion),
                     sum(edge_valence[e] for e in adj_boson))
                    for adj_fermion, adj_boson in adjacence)
        if any(val not in allowed_valencies for val in valences):
            continue

        # check that legs are correct
        if ext_fermion:
            fermion_legs = sum(1 for adj in ext_adj for e in adj
                               if full_weights[e] == fermion)
            if fermion_legs != ext_fermion:
                continue
        if ext_boson:
            boson_legs = sum(1 for adj in ext_adj for e in adj
                             if full_weights[e] == boson)
            if boson_legs != ext_boson:
                continue

        final_weights = tuple(boson if e in boson_edges else fermion
                              for e in graph.edges_set)

        # handling the fermion directions
        fermion_adj = [a for a, _ in adjacence]

        for fermion_weights in product((-1, 1), repeat=len(fermion_edges)):

            weight_table = {}
            pos = 0
            for i in graph.edges_set:
                if i not in boson_edges:
                    weight_table[i] = fermion_weights[pos]
                    pos += 1

            # fermion edges must form oriented cycles
            if any(sum(dir_sign(e, v) * weight_table[e]
                       for e in adj if not is_sl[e])
                   for v, adj in zip(int_vtcs, fermion_adj) if adj):
                continue

            edges = []
            pos = 0
            for i, edge in enumerate(graph.edges):
                if i in boson_edges:
                    edges.append(edge)
                else:
                    edges.append(edge if fermion_weights[pos] == 1
                                 else flip(edge))
                    pos += 1

            g = WeightedGraph(edges, final_weights)

            # what happens below for external paths ?
            _, cycles = g.cycle_decomposition(fermion_edges)
            if any(len(c) % 2 for c in cycles):
                continue
            yield g
