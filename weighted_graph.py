
"""weighted_graph.py: This file is part of the feyncop/feyngen package.
    Implements the WeightedGraph class. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

from math import factorial
import copy
from itertools import permutations

from stuff import double_factorial, flip
from graph import Graph


wDict = ['0', 'f', 'A', 'c']
# what is the convention here ?
# ?, fermion, boson, ?


class WeightedGraph(Graph):
    """This class extends the basic utilities in the Graph class by the tools
        to handle QED and Yang-Mills graphs."""

    def __init__(self, edges, edge_weights=None, symmetry_factor=0):
        """
        Initialize the WeightedGraph class.

        INPUT:

        - either a weighted graph

        - or (edges, edge_weights)

        - or (edges, edge_weights, symmetry factor)

        EXAMPLES::

            sage: G = WeightedGraph([[0,1],[0,1],[0,1],[0,1]],[0,1,2,3]); G
            G[[0,1,0],[0,1,f],[0,1,A],[0,1,c]]

            sage: WeightedGraph([[0,1],[1,2]],[0,1],2)
            G[[0,1,0],[1,2,f]]/2
        """
        if isinstance(edges, WeightedGraph):
            # input = (graph,)
            G = edges
            edges = G.edges
            edge_weights = G.edge_weights
            symmetry_factor = G.symmetry_factor
        else:
            # input = (edges, edge_weights, maybe integer)
            if edge_weights is None:
                raise TypeError("weights are missing")
            if len(edges) != len(edge_weights):
                raise ValueError("bad number of weights")

        super().__init__(edges, symmetry_factor)
        self.edge_weights = edge_weights

    def get_edge_str(self, e):
        """Return a readable string of the edges of the graph."""

        v1, v2 = self.edges[e]
        w = self.edge_weights[e]
        return "[%d,%d,%c]" % (v1, v2, wDict[w])

    def sage(self):
        """
        Transform into a Sage graph.

        EXAMPLES::

            sage: from weighted_graph import *
            sage: WeightedGraph([[0,1],[1,2]],[1,2]).sage()
            Looped multi-graph on 3 vertices
        """
        from sage.graphs.graph import Graph as SageGraph
        edges = [(a, b, c) for (a, b), c in zip(self.edges, self.edge_weights)]
        return SageGraph(edges, loops=True, multiedges=True,
                         format="list_of_edges")

    def get_edges_tuple(self):
        """Get a unique tuple to identify the graph. (Unique only for every labeling)."""

        return tuple(sorted((tuple(sorted(edge) if w == 2 else edge), w)
                            for edge, w in zip(self.edges, self.edge_weights)))

    def graph_from_sub_edges(self, sub_edges):
        """Create a new graph from a sub set of its edges."""

        sub_graph = super().graph_from_sub_edges(sub_edges)
        sub_graph.edge_weights = tuple(self.edge_weights[e]
                                       for e in sorted(sub_edges))
        return sub_graph

    def delete_one_vertex(self):
        """
        Create a weighted graph with one arbitrary vertex removed.

        EXAMPLES::

            sage: G = WeightedGraph([[0,1],[0,1],[0,1],[0,2],[1,3],[2,3],[2,3]],[1,1,1,1,1,2,2])
            sage: G.delete_one_vertex()
            G[[0,1,f],[0,1,f],[0,1,f],[0,2,f]]
        """
        edges = self.edges
        weights = self.edge_weights
        n = self.num_verts - 1
        data = [(e, w) for e, w in zip(edges, weights)
                if n not in e]
        return WeightedGraph([p[0] for p in data], [p[1] for p in data])

    def sub_edges_by_weight(self, weight):
        """Returns all subedges with a certain weight."""

        return frozenset(e for e, w in enumerate(self.edge_weights)
                         if w == weight)

    @property
    def residue_type(self):
        """Returns the residue type of the graph."""

        def dir_e(e, v):
            if self.edge_weights[e] == 2:
                return 1
            return -1 if v == self.edges[e][0] else 1

        ext_types = (dir_e(e, v) * self.edge_weights[e]
                     for v in self.external_vtcs_set
                     for e in self.adj_edges(v, self.edges_set))
        return tuple(sorted(ext_types))

    def get_vtx_type(self, v):
        """Returns the type of the vertex v in the same format as
            residue_type."""

        def dir1(e, v):
            if self.edge_weights[e] == 2:
                return 1
            return -1 if v == self.edges[e][0] else 1

        def dir2(e, v):
            if self.edge_weights[e] == 2:
                return 1
            return 1 if v == self.edges[e][0] else -1

        adj_types = [dir1(e, v) * self.edge_weights[e]
                     for e in self.adj_edges(v, self.edges_set)]
        adj_types.extend(dir2(e, v) * self.edge_weights[e]
                         for e in self.edges_set if self.edges[e] == (v, v))

        return tuple(sorted(adj_types))

    def get_vtcs_coloring(self):
        """Helper function: Calculate the vertex coloring in a format suitable
            for the canonical labeling calculation."""

        # All vertices with different numbers of selfloops of different type
        # are colored in another way.
        dictWeights = {edge: self.edge_weights[e] for e, edge in enumerate(self.edges)}
        edge_degree_counter = self.edge_degree_counter(self.edges_set)
        selfloop_degree_list = [(edge_degree_counter[(v, v)], dictWeights[(v, v)]
                                 if edge_degree_counter[(v, v)] else 2)
                                for v in self.internal_vtcs_set]

        # Sorting is important for the v even for all similar mul!
        selfloop_multiplicity_list = sorted((mul, v)
                                            for v, mul in zip(self.internal_vtcs_set, selfloop_degree_list))
        max_selfloop_multiplicity, _ = selfloop_multiplicity_list[-1][0] if selfloop_multiplicity_list else (0, 2)

        self_loop_list = [frozenset(vtx
                                    for (mul, we), vtx in selfloop_multiplicity_list
                                    if mul == i and we == w)
                          for i in range(max_selfloop_multiplicity + 1)
                          for w in (1, 2, 3)]

        # External vertices all have the same color still.
        return self_loop_list + [self.external_vtcs_set]

    def get_edges_coloring(self, edges_set):
        """Helper function: Calculate the edge coloring in a format suitable
            for the canonical labeling calculation."""

        # Fermions, bosons and ghosts need different color classes.
        fermion_edges_set = self.sub_edges_by_weight(1) & edges_set
        boson_edges_set = self.sub_edges_by_weight(2) & edges_set
        ghost_edges_set = self.sub_edges_by_weight(3) & edges_set

        fermion_edges = frozenset(self.edges[i] for i in fermion_edges_set if not self.is_selfloop(self.edges[i]))
        ghost_edges = frozenset(self.edges[i] for i in ghost_edges_set if not self.is_selfloop(self.edges[i]))
        boson_edges = frozenset(self.edges[i] for i in boson_edges_set)

        # Fermions and ghosts need orientation. Bosons not!
        # For higher performance some special cases of boson-fermion-ghost
        # edge combinations are included.
        def normalize(edge):
            return (max(edge), min(edge))

        fermion_loops = frozenset(normalize(edge) for edge in fermion_edges if flip(edge) in fermion_edges)
        ghost_loops = frozenset(normalize(edge) for edge in ghost_edges if flip(edge) in ghost_edges)
        reduced_fermion_edges = fermion_edges - fermion_loops - frozenset(flip(edge) for edge in fermion_loops)
        reduced_ghost_edges = ghost_edges - ghost_loops - frozenset(flip(edge) for edge in ghost_loops)
        boson_fermion_loops = frozenset(edge for edge in reduced_fermion_edges if flip(edge) in boson_edges or edge in boson_edges)
        boson_ghost_loops = frozenset(edge for edge in reduced_ghost_edges if flip(edge) in boson_edges or edge in boson_edges)

        reduced_boson_edges = boson_edges - boson_fermion_loops - frozenset(flip(edge) for edge in boson_fermion_loops) - boson_ghost_loops - frozenset(flip(edge) for edge in boson_ghost_loops)

        dbl_boson_edges = reduced_boson_edges | frozenset(flip(edge) for edge in reduced_boson_edges)

        if dbl_boson_edges & reduced_fermion_edges or \
           dbl_boson_edges & reduced_ghost_edges:
            print(dbl_boson_edges, reduced_fermion_edges)
            raise RuntimeError

        # Calculate the boson coloring as in the Graph class.
        boson_coloring = super().get_edges_coloring(boson_edges_set)

        return [dbl_boson_edges | reduced_fermion_edges | reduced_ghost_edges,
                fermion_loops, boson_fermion_loops,
                ghost_loops, boson_ghost_loops,
                reduced_ghost_edges - boson_ghost_loops] + boson_coloring[1:]

    def get_trivial_symmetry_factor(self):
        """Calculates the trivial factor in the symmetry factor. Only
            considers edge multiplicity and self loops."""

        grpSize = 1
        boson_edges = self.sub_edges_by_weight(2)
        edge_degree_counter = self.edge_degree_counter(boson_edges)
        for mul_edge_deg in (m for edge, m in edge_degree_counter.items() if not self.is_selfloop(edge)):
            grpSize *= factorial(mul_edge_deg)

        for selfloop_deg in (m for edge, m in edge_degree_counter.items() if self.is_selfloop(edge)):
            grpSize *= double_factorial(2 * selfloop_deg)
        return grpSize

    def permute_external_edges(self):
        """Generate all possible graphs with fixed external legs from the
            graph provided that the graph is non-leg-fixed."""
        class FixedGraph(type(self)):
            def get_vtcs_coloring(self):
                vtcs_coloring = super().get_vtcs_coloring()

                vtcs_coloring = [c - self.external_vtcs_set for c in vtcs_coloring]
                vtcs_coloring.extend(frozenset([v]) for v in sorted(self.external_vtcs_set))

                return vtcs_coloring

        extern_boson_vtcs = \
            frozenset(v for e in self.sub_edges_by_weight(2) for v in self.edges[e]) \
            & self.external_vtcs_set
        extern_in_fermion_vtcs = \
            frozenset(self.edges[e][0] for e in self.sub_edges_by_weight(1)) \
            & self.external_vtcs_set
        extern_out_fermion_vtcs = \
            frozenset(self.edges[e][1] for e in self.sub_edges_by_weight(1)) \
            & self.external_vtcs_set
        extern_in_ghost_vtcs = \
            frozenset(self.edges[e][0] for e in self.sub_edges_by_weight(3)) \
            & self.external_vtcs_set
        extern_out_ghost_vtcs = \
            frozenset(self.edges[e][1] for e in self.sub_edges_by_weight(3)) \
            & self.external_vtcs_set

        extern_vtcs_list = list(extern_boson_vtcs) + \
            list(extern_in_fermion_vtcs) + \
            list(extern_out_fermion_vtcs) + \
            list(extern_in_ghost_vtcs) + \
            list(extern_out_ghost_vtcs)
        if frozenset(extern_vtcs_list) != self.external_vtcs_set:
            raise

        vtcs_list = list(self.internal_vtcs_set) + \
            extern_vtcs_list

        def relabel_edge(v12, m):
            v1, v2 = v12
            return (m[v1], m[v2])

        for perm0 in permutations(extern_boson_vtcs):
            for perm1 in permutations(extern_in_fermion_vtcs):
                for perm2 in permutations(extern_out_fermion_vtcs):
                    for perm3 in permutations(extern_in_ghost_vtcs):
                        for perm4 in permutations(extern_out_ghost_vtcs):

                            new_vtcs_list = tuple(self.internal_vtcs_set) + \
                                perm0 + perm1 + perm2 + perm3 + perm4
                            m = dict(zip(vtcs_list, new_vtcs_list))

                            yield FixedGraph(
                                [relabel_edge(edge, m) for edge in self.edges], self.edge_weights, 0)

    @property
    def clean_graph(self):
        """Orders the edge- and weight list of the graph in a transparent manner."""

        ext_sorter = (e in self.external_edges_set for e, edge in enumerate(self.edges))

        def norm(edge):
            return (max(edge), min(edge))

        edges = [norm(edge) if w == 2 else edge
                 for w, edge in zip(self.edge_weights, self.edges)]
        xwe_list = sorted(zip(ext_sorter, self.edge_weights, edges))
        edges = [edge for _, _, edge in xwe_list]
        weights = [w for _, w, _ in xwe_list]
        g = copy.copy(self)
        g.edges = tuple(edges)
        g.edge_weights = tuple(weights)
        g.prepare_graph()

        return g
