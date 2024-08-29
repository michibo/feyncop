
""" graph.py: This file is part of the feyncop/feyngen package.
    Implements the Graph class. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

from math import *
from stuff import *
import collections
import itertools
import copy
import nauty_wrapper


class Graph:
    """This class incorporates all the graph theoretic tools necessary for basic
        Feynman graph handling."""

    def __init__(self, edges, symmetry_factor=0):
        """Initializes the Graph class. Edges and symmetry_factor can be provided."""

        self.edges = edges
        self.symmetry_factor = symmetry_factor

        self.prepare_graph()

    def __str__(self):
        g_string = ",".join(self.get_edge_str(e) for e in self.edges_set)
        if self.symmetry_factor:
            return "G[%s]/%d" % (g_string, self.symmetry_factor)
        else:
            return "G[%s]" % (g_string)
    __repr__ = __str__

    def __eq__(self, other):
        """Compare two graphs. Only the labelings are compared. To check
            non-isomorphy both graphs must be canonically labeled."""

        return self.get_edges_tuple() == other.get_edges_tuple()

    def __ne__(self, other):
        return not self == other

    def __hash__(self):
        """Create hash from the labeling. The hash is unique for the
            isomorphism class of the graph if the labeling is
            canonical."""

        return hash(self.get_edges_tuple())

    def prepare_graph(self):
        """Initializes dictionaries and sets for high performance graph handling."""

        self.edges_set = frozenset(range(len(self.edges)))
        vtcs = frozenset(v for edge in self.edges for v in edge)
        self.valency_dict = collections.defaultdict(lambda: 0)
        for v1,v2 in self.edges:
            self.valency_dict[v1] += 1
            self.valency_dict[v2] += 1

        self.internal_vtcs_set = frozenset(v for v in vtcs if self.valency_dict[v] > 1)
        self.external_vtcs_set = frozenset(v for v in vtcs if self.valency_dict[v] == 1)
        self.internal_edges_set = frozenset(e for e,(v1,v2) in enumerate(self.edges) if v1 in self.internal_vtcs_set and v2 in self.internal_vtcs_set)
        self.external_edges_set = self.edges_set - self.internal_edges_set

    def get_edge_str(self, e):
        """Return a readable string of the edges of the graph."""

        v1, v2 = self.edges[e]
        return "[%d,%d]" % (v1, v2)

    def get_edges_tuple(self):
        """Get a unique tuple to identify the graph. (Unique only for every labeling)."""

        return tuple(sorted(tuple(sorted(edge)) for edge in self.edges))

    def graph_from_sub_edges(self, sub_edges):
        """Create a new graph from a sub set of its edges."""

        edges = tuple(self.edges[e] for e in sorted(sub_edges))

        sub_graph = copy.copy(self)
        sub_graph.edges = edges
        sub_graph.prepare_graph()
        sub_graph.symmetry_factor = 0

        return sub_graph

    def adj_edges(self, v, sub_edges):
        """Return the adjacent edges to vertex v. Only consider edges in
            sub_edges"""

        is_adj = lambda xy: (xy[0] == v) or (xy[1] == v)
        return (e for e, edge in enumerate(self.edges)
                if e in sub_edges and is_adj(edge))

    def vtx_valence(self, v, sub_edges):
        """Return the valence of vertex v. Only consider edges in sub_edges as
            relevant."""

        return sum(1 for e in sub_edges for w in self.edges[e] if v == w)

    def is_selfloop(self, edge):
        """Is the edge a selfloop?"""

        x,y = edge
        return x == y

    def vtcs_set_sub_edges(self, sub_edges):
        """Return all the vertices adjacent to the edges in
            sub_edges."""

        return frozenset(v for e in sub_edges for v in self.edges[e])

    def edge_degree_counter(self, sub_edges):
        """Return a counter of edge multiplicity."""

        def norm(xy):
            x, y = xy
            return xy if x > y else (y, x)

        return collections.Counter(norm(self.edges[e]) for e in sub_edges)

    def dfs(self, v, sub_edges, back_edge_visitor, forward_edge_visitor, discovered=set(), forward_edges=set(), back_edges=set(), trace=(), trace_edges=()):
        """Implements the dfs algorithm on the graph in a simple recursive manner.
            v is the start vertex. Only sub_edges are considered. The other
            variables carry the standard names. See in Corman, Leiserson,
            Rivest, Stein - Algorithms for details. """

        discovered|= set([v])
        adj_edges = frozenset(self.adj_edges(v, sub_edges))
        get_adj = lambda xy: xy[1] if xy[0] == v else xy[0]

        for i in adj_edges:
            if i in forward_edges or i in back_edges:
                continue

            adj_v = get_adj(self.edges[i])
            if adj_v in discovered:
                back_edges|= set([i])
                if back_edge_visitor:
                    back_edge_visitor(v, adj_v, i, trace, trace_edges)
            else:
                forward_edges|= set([i])
                if forward_edge_visitor:
                    forward_edge_visitor(v, adj_v, i, trace, trace_edges)

                self.dfs(adj_v, sub_edges, back_edge_visitor, forward_edge_visitor, discovered, forward_edges, back_edges, trace + (v,), trace_edges + (i,))

    def cycle_decomposition(self, sub_edges):
        """Decomposes the graph into its cycles using the dfs function. Only
            sub_edges are considered."""

        cycles, cycles_edges = [], []

        def back_edge_visitor(v, adj_v, back_edge, trace, trace_edges):
            if adj_v == v:
                cycles.append((v,))
                cycles_edges.append((back_edge,))
            else:
                index = next(j for j,vtx in enumerate(trace) if vtx == adj_v)
                cycles.append(trace[index:] + (v,))
                cycles_edges.append(trace_edges[index:] + (back_edge,))

        def forward_edge_visitor(v, adj_v, forward_edge, trace, trace_edges):
            pass

        sub_edges_cpy = set(sub_edges)
        while sub_edges_cpy:
            discovered, f_edges, b_edges = set(), set(), set()
            v,_ = self.edges[next(iter(sub_edges_cpy))]
            self.dfs(v, sub_edges_cpy, back_edge_visitor, forward_edge_visitor, discovered, f_edges, b_edges)
            sub_edges_cpy -= f_edges | b_edges

        return cycles, cycles_edges

    def is_connected_sub_edges(self, sub_edges):
        """Returns true if the sub graph consisting of sub_edges is
            connected."""

        if not sub_edges:
            if len(self.internal_vtcs_set) <= 1:
                return True
            return False

        vtcs_set = self.vtcs_set_sub_edges(sub_edges)

        f_edges, b_edges = (set(), set())

        try:
            v,_ = self.edges[next(iter(sub_edges))]
        except StopIteration:
            print("Error: Graph is not connected?")
            raise

        # Connectedness is checked using the dfs algorithm.
        discovered = set()
        self.dfs(v, sub_edges, None, None, discovered, f_edges, b_edges)
        return sub_edges - (f_edges | b_edges) == set() and vtcs_set - discovered == set()

    @property
    def is_connected(self):
        """True if the whole graph is connected."""
        return self.is_connected_sub_edges(self.edges_set)

    def cntd_components_sub_edges(self, sub_edges):
        """Generates the connected components of the subgraph
            consisting of sub_edges."""

        pre_discovered = set()

        left_edges = set(sub_edges)
        while left_edges:
            v,_ = self.edges[next(iter(left_edges))]
            discovered, f_edges, b_edges = set(), set(), set()

            # Use dfs iteratedly until all "islands" are identified.
            self.dfs(v, sub_edges, None, None, discovered, f_edges, b_edges)

            if discovered & pre_discovered:
                print("Warning: Connected components error")
                raise
            pre_discovered |= discovered

            component_edges = f_edges | b_edges
            left_edges -= component_edges

            yield component_edges

    @property
    def cntd_components(self):
        """Generates the connected components of the graph."""

        yield from self.cntd_components_sub_edges(self.edges_set)

    def is_edge_2_connected_sub_edges(self, sub_edges):
        """True if the subgraph consisting of sub_edges is edge-2-connected."""

        if not self.is_connected_sub_edges(sub_edges):
            return False

        # Use cycle decomposition to check that no bridge-edges are in the
        # graph.
        cycles, cycles_edges = self.cycle_decomposition(sub_edges&self.internal_edges_set)
        loop_edges = set(e for cycle in cycles_edges for e in cycle)
        nonloop_edges = sub_edges&self.internal_edges_set - loop_edges

        return nonloop_edges == set() and loop_edges

    @property
    def is_edge_2_connected(self):
        """True if the graph is edge-2-connected."""

        return self.is_edge_2_connected_sub_edges(self.edges_set)

    @property
    def is_vtx_2_connected(self):
        """True if the graph is vertex-2-connected."""

        if not self.is_connected:
            return False

        # Use cycle decomposition to check that no cut-vertices are in the
        # graph.
        cycles, cycles_edges = self.cycle_decomposition(self.internal_edges_set)

        bicntd_components = [set(cycle) for cycle in cycles_edges]
        for i in range(len(bicntd_components)-1):
            c = bicntd_components[i]
            bicntd_components[i+1:] = [c|d if c&d else d for d in bicntd_components[i+1:]]

        return bicntd_components and bicntd_components[-1] == self.internal_edges_set and len(bicntd_components[-1]) > 1

    @property
    def is_tadpole(self):
        """True if the graph is a tadpole."""

        # For not being a tadpole grpah there must not be any selfloops and
        # all the biconnected components must be connected to an external
        # vertex."""

        edge_degree_counter = self.edge_degree_counter(self.edges_set)
        selfloop_degree_list = [edge_degree_counter[(v,v)] for v in self.internal_vtcs_set]
        if any(tp_deg>0 for tp_deg in selfloop_degree_list):
            return True

        if len(self.external_vtcs_set) == 1:
            return True

        cps0_len = sum(1 for comp in self.cntd_components)
        for v in self.internal_vtcs_set:
            sub_edges = self.edges_set - set(self.adj_edges(v, self.edges_set))
            components = list(self.cntd_components_sub_edges(sub_edges))

            for comp in components:
                if all(e not in self.external_edges_set for e in comp):
                    return True

        return False

    def get_vtcs_coloring(self):
        """Helper function: Calculate the vertex coloring in a format suitable
            for the canonical labeling calculation."""

        edge_degree_counter = self.edge_degree_counter(self.edges_set)
        # All vertices with different numbers of selfloops are colored in
        # another way.
        selfloop_degree_list = [edge_degree_counter[(v,v)] for v in self.internal_vtcs_set]

        # Sorting is important for the v even for all similar mul!
        selfloop_multiplicity_list = sorted((mul,v) for v, mul in zip(self.internal_vtcs_set, selfloop_degree_list))
        (max_selfloop_multiplicity, _) = selfloop_multiplicity_list[-1] if selfloop_multiplicity_list else (0, 0)

        self_loop_list = [frozenset(vtx
                                    for mul, vtx in selfloop_multiplicity_list
                                    if mul == i)
                          for i in range(max_selfloop_multiplicity + 1)]

        # External vertices all have the same color still.
        return self_loop_list + [self.external_vtcs_set]

    def get_edges_coloring(self, edges_set):
        """Helper function: Calculate the edge coloring in a format suitable
            for the canonical labeling calculation."""

        # Selfloops don't need color.
        non_selfloop_edges_set = frozenset(e for e in self.edges_set if not self.is_selfloop(self.edges[e])) & edges_set

        edge_degree_counter = self.edge_degree_counter(non_selfloop_edges_set)

        # Edges with different multiplicity get different colors.
        edge_multiplicity_list = sorted((mul, edge) for edge, mul in edge_degree_counter.items())
        (max_edge_multiplicity, _) = edge_multiplicity_list[-1] if edge_multiplicity_list else (0, 0)

        edge_coloring = [[edge for mul, edge in edge_multiplicity_list
                         if mul == i]
                         for i in range(1, max_edge_multiplicity + 1)]

        if edge_coloring:
            def flip(xy):
                x, y = xy
                return (y, x)
            edge_coloring[0] += [flip(edge) for edge in edge_coloring[0]]

        return edge_coloring

    def nygraph(self, colored_vtcs, colored_edges):
        """Helper function: Calculate the simple graph without selfloops
            that can be used as input for nauty. The edge and vertex
            colorings are transformed into a single vertex coloring
            suitable as input for nauty."""

        vp_list = [v for part in colored_vtcs for v in part]
        ny_lab = dict((v,i) for i,v in enumerate(vp_list))

        colored_vtcs = [[ny_lab[v] for v in part] for part in colored_vtcs]
        colored_edges = [[(ny_lab[v1],ny_lab[v2]) for v1,v2 in edgeset] for edgeset in colored_edges]

        ny_edges = colored_edges[0] if colored_edges else []

        # Transform edge coloring to vertex coloring by adding extra vertices.
        num_vertices = len(vp_list)
        for edge_set in colored_edges[1:]:
            if not edge_set:
                continue

            new_vertices = frozenset(range(num_vertices, num_vertices + len(edge_set)))
            new_edges = [(src, help_vtx) for help_vtx, (src, tgt) in zip(new_vertices, edge_set)] + \
                        [(tgt, help_vtx) for help_vtx, (src, tgt) in zip(new_vertices, edge_set)]

            colored_vtcs.append(new_vertices)
            ny_edges.extend(new_edges)

            num_vertices += len(new_vertices)

        return (colored_vtcs, ny_edges, num_vertices, ny_lab)

    def get_canonical_edges(self, colored_vtcs, colored_edges):
        """Calculates canonically labeled edges and the size of
            the automorphism group. """

        if not self.edges:
            return tuple(self.edges), 1

        # Calc the simple graph.
        (ny_colored_vtcs, ny_edges, num_ny_vertices, ny_lab) = self.nygraph(colored_vtcs, colored_edges)

        # Give it to nauty.
        (lab, grpSize, orbits) = nauty_wrapper.get_canonical_labeling(num_ny_vertices, ny_edges, ny_colored_vtcs)

        vtx_list = [v for part in ny_colored_vtcs for v in part]
        m = dict((old, new) for old,new in zip(lab, vtx_list))

        def relabel_edge(v12):
            v1, v2 = v12
            return (m[ny_lab[v1]], m[ny_lab[v2]])

        # Relabel the edges according to the canonical labeling.
        canonical_edges = tuple(relabel_edge(edge) for edge in self.edges)
        return canonical_edges, grpSize

    def get_trivial_symmetry_factor(self):
        """Calculates the trivial factor in the symmetry factor. Only
            considers edge multiplicity and self loops."""

        grpSize = 1
        edge_degree_counter = self.edge_degree_counter(self.edges_set)
        for mul_edge_deg in (m for edge, m in edge_degree_counter.items() if not self.is_selfloop(edge)):
            grpSize*= factorial(mul_edge_deg)

        for selfloop_deg in (m for edge, m in edge_degree_counter.items() if self.is_selfloop(edge)):
            grpSize*= double_factorial(2*selfloop_deg)
        return grpSize

    @property
    def unlabeled_graph(self):
        """Returns an unlabeled/canonically labeled graph."""
        vtcs_coloring = self.get_vtcs_coloring()
        edges_coloring = self.get_edges_coloring(self.edges_set)

        canonical_edges, grpSize = self.get_canonical_edges(vtcs_coloring, edges_coloring)
        grpSize*= self.get_trivial_symmetry_factor()

        unlabeled_graph = copy.copy(self)
        unlabeled_graph.edges = canonical_edges
        unlabeled_graph.symmetry_factor = grpSize

        return unlabeled_graph.clean_graph

    def permute_external_edges(self):
        """Generate all possible graphs with fixed external legs from the
            graph provided that the graph is non-leg-fixed."""

        # Fixed graphs get their own child class, because the
        # canonical labeling calculation must be altered.
        class FixedGraph(type(self)):
            def get_vtcs_coloring(self):
                vtcs_coloring = super().get_vtcs_coloring()

                vtcs_coloring = [c - self.external_vtcs_set for c in vtcs_coloring]
                vtcs_coloring.extend(frozenset([v]) for v in sorted(self.external_vtcs_set))

                return vtcs_coloring

        for perm in itertools.permutations(sorted(self.external_vtcs_set)):
            m = tuple(self.internal_vtcs_set) + perm

            def relabel_edge(v12):
                v1, v2 = v12
                return (m[v1], m[v2])

            yield FixedGraph(
                [relabel_edge(edge) for edge in self.edges], 0)

    @property
    def clean_graph(self):
        """Orders the edge list of the graph in a transparent manner."""

        ext_sorter = (e in self.external_edges_set for e,edge in enumerate(self.edges))

        norm = lambda edge: (max(edge), min(edge))
        edges = [norm(edge) for edge in self.edges]
        xe_list = sorted(zip(ext_sorter, edges))
        edges = [edge for x,edge in xe_list]
        g = copy.copy(self)
        g.edges = tuple(edges)

        g.prepare_graph()

        return g
