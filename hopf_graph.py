
"""hopf_graph.py: This file is part of the feyncop/feyngen package.
    Implements the HopfGraph class. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

import collections
import copy
import itertools
import sys

from graph import Graph
from weighted_graph import WeightedGraph


class HopfGraph(WeightedGraph):
    """This class incorporates all the tools necessary for
        the treatment of graphs in the scope of the
        Hopf algebra of Feynman graphs."""

    def get_graph_str(self, ym=False):
        """Get the string representation of the graph either as
            Yang-Mills/QED or phi^k graph."""

        if ym:
            return str(self)
        return str(Graph(self.edges, self.symmetry_factor))

    def is_primitive(self, dimension, ym=False):
        """True if the graph is primitive."""

        try:
            _ = next(self.reduced_coproduct(dimension, ym))
        except StopIteration:
            # True if the reduced coproduct is 0.
            return True

        return False

    def reduced_coproduct(self, dimension, ym=False):
        """Calculate the reduced coproduct of the graph in the given
            dimension."""

        if not self.is_edge_2_connected:
            print("Warning: Calculation of non 1PI graph coproduct omitted:", self, file=sys.stderr)
            return

        yield from self.reduced_coproduct_on_edges(dimension, self.internal_edges_set, ym)

    def reduced_coproduct_with_residue_graph(self, dimension, ym=False):
        """Calculate the reduced coproduct of the graph in the given
            dimension and the corresponding cograph."""

        for sub_edges_list in self.reduced_coproduct(dimension, ym):
            sub_edges_set = frozenset(e for sub_edges in sub_edges_list for e in sub_edges)
            residue_edges = self.calc_residue_edges(self.edges_set, sub_edges_set)

            residue_graph = copy.copy(self)
            residue_graph.edges = tuple(residue_edges)
            residue_graph.prepare_graph()

            yield sub_edges_list, residue_graph

    def eval_subedges_for_reduced_coproduct(self, sub_edges, dimension, sub_edges_list, ym=False):
        """Helper function: Evaluates omega_D on a given subgraph and checks if it is
            superficially divergent."""

        for component in self.cntd_components_sub_edges(sub_edges):
            num_vtcs = len(self.vtcs_set_sub_edges(component))

            num_loops = len(component) - num_vtcs + 1
            edge_weights = (self.edge_weights[e] for e in component)
            denom = sum(w if w != 3 else 2 for w in edge_weights)
            if ym:
                vtx_types = (self.get_vtx_type(v) for v in self.vtcs_set_sub_edges(component))
                # The ghost vertices and the triple gluon vertices get an extra -1.
                denom -= sum(1 for t in vtx_types if (t == (2, 2, 2)) or (t == (-3, 2, 3)))

            dim_int = dimension * num_loops
            omega = denom - dim_int

            if omega > 0:
                return False

            if not self.is_edge_2_connected_sub_edges(component):
                return False

            sub_edges_list.append(component)

        return True

    def reduced_coproduct_on_edges(self, dimension, edges_set, ym=False):
        """Calculate the reduced coproduct in the given dimension
            only including edgers in edges_set."""

        for i in range(1, len(edges_set)):
            for sub_edges in itertools.combinations(edges_set, i):
                sub_edges_list = []
                if self.eval_subedges_for_reduced_coproduct(sub_edges, dimension, sub_edges_list, ym):
                    yield sub_edges_list

    def reduced_coproduct_unlabeled(self, dimension, ym=False):
        """Calculate the reduced coproduct and give the result
            as tensor products of unlabeled graphs."""

        for sub_edges_list, residue_graph in self.reduced_coproduct_with_residue_graph(dimension, ym):
            sub_graphs = collections.Counter()
            for sub_edges in sub_edges_list:
                if not sub_edges:
                    continue
                sub_graph = self.sub_graph_with_legs(sub_edges)
                sub_graphs[sub_graph.unlabeled_graph] += 1

            sub_edges_set = frozenset(e for sub_edges in sub_edges_list for e in sub_edges)
            residue_edges_set = (self.edges_set - sub_edges_set) | self.external_edges_set
            residue_graph_reduced = residue_graph.graph_from_sub_edges(residue_edges_set)
            residue_graph_reduced.clean_of_val2_vtcs()

            yield sub_graphs, residue_graph_reduced.unlabeled_graph

    def calc_residue_edges(self, edges_set, sub_edges):
        """Helper function: Calculates the residue edges corresponding to a
            certain subgraph."""

        residue_edges = list(self.edges)

        m = {v: v for v in frozenset(v for edge in self.edges for v in edge
                                     if v != -1)}
        m[-1] = -1
        for sub_edge in sub_edges:
            w1, w2 = self.edges[sub_edge]
            r = m[w1], m[w2]
            o, n = max(r), min(r)

            m = {v: (n if m[v] == o else m[v]) for v in m}

            residue_edges[sub_edge] = (-1, -1)
            residue_edges = [(m[v1], m[v2]) for v1, v2 in residue_edges]
        return residue_edges

    def sub_graph_with_legs(self, sub_edges):
        """Helper function: Adds legs to a given subgraph in a canonical way."""

        not_edges = self.edges_set - sub_edges
        vtcs = self.vtcs_set_sub_edges(sub_edges)
        not_vtcs = self.vtcs_set_sub_edges(not_edges)
        ext_vtcs = vtcs & not_vtcs

        def adj_vtx_e(e, v):
            x, y = self.edges[e]
            return y if x == v else x

        ext_edges = [(adj_vtx_e(e, v), v, e) for v in ext_vtcs
                     for e in self.adj_edges(v, not_edges)]
        ext_edges += [(v, v, e) for v in ext_vtcs for e in self.adj_edges(v, not_edges)
                      if adj_vtx_e(e, v) == v]
        ext_vtcs_offset = max(vtcs) + 1

        def dir_e(vn, v, e):
            if self.edges[e][0] == v:
                return (v, vn)
            return (vn, v)

        new_edges = [(dir_e(i + ext_vtcs_offset, v, e), e)
                     for i, (vx, v, e) in enumerate(sorted(ext_edges))]
        sub_graph = self.graph_from_sub_edges(sub_edges)
        sub_graph.edges = tuple(list(sub_graph.edges) + [edge for edge, e in new_edges])
        sub_graph.edge_weights = tuple(list(sub_graph.edge_weights) + [self.edge_weights[e] for edge, e in new_edges])
        sub_graph.prepare_graph()

        return sub_graph

    def clean_of_val2_vtcs(self):
        """Helper function: Removes the valency 2 vertices in the graph."""

        while True:
            bad_vtcs = [v for v in self.internal_vtcs_set if self.vtx_valence(v, self.edges_set) == 2 and (v, v) not in self.edges]
            if not bad_vtcs:
                break
            bad_vtcs = [bad_vtcs[0]]
            pre_bad_edges = {v: tuple(self.adj_edges(v, self.edges_set)) for v in bad_vtcs}
            bad_edges = frozenset(e for v, edges in pre_bad_edges.items() for e in edges)

            def gen_new_edges():
                for v, edges in pre_bad_edges.items():
                    def adj_v(xy):
                        return xy[0] if xy[1] == v else xy[1]
                    e1, e2 = edges
                    w1, w2 = self.edge_weights[e1], self.edge_weights[e2]
                    if w1 != w2:
                        raise
                    v1, v2 = adj_v(self.edges[e1]), adj_v(self.edges[e2])
                    if v1 != self.edges[e1][0]:
                        v1, v2 = v2, v1

                    yield (v1, v2), w1

            good_edges_w = list(gen_new_edges())

            new_edges = [self.edges[e] for e in self.edges_set - bad_edges] + [edge for edge, w in good_edges_w]
            new_weights = [self.edge_weights[e] for e in self.edges_set - bad_edges] + [w for edge, w in good_edges_w]
            self.edges = tuple(new_edges)
            self.edge_weights = tuple(new_weights)
            self.prepare_graph()
