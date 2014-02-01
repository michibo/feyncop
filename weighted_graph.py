
from math import *
import copy, itertools
import collections

from stuff import *
from graph import Graph

class WeightedGraph(Graph):
    def __init__( self, edges, edge_weights, symmetry_factor=0 ):
        if len(edges) != len(edge_weights):
            raise
        
        super(WeightedGraph, self).__init__( edges, symmetry_factor )
        self.edge_weights = edge_weights

    def sub_edges_by_weight( self, weight ):
        return frozenset( e for e,w in enumerate(self.edge_weights) if w == weight )

    def graph_from_sub_edges( self, sub_edges ):
        sub_graph = super(WeightedGraph, self).graph_from_sub_edges( sub_edges )
        sub_graph.edge_weights = tuple( self.edge_weights[e] for e in sorted(sub_edges) )

        return sub_graph

    def get_edges_tuple( self ):
        return tuple( sorted( ( tuple( sorted(edge) if w!=1 else edge ), w) for edge,w in zip(self.edges,self.edge_weights) ) )

    def get_edge_str( self, e ):
        #if all( w == 2 for w in self.edge_weights ):
        #    return super(WeightedGraph, self).get_edge_str(e)

        v1,v2 = self.edges[e]
        w = self.edge_weights[e]
        wDict = [ '0', 'f', 'A', 'c' ]
        return "[%d,%d,%c]" % (v1,v2,wDict[w])

    def get_graph_sign( self ):
        fermion_edges_set = self.sub_edges_by_weight(1)
        cycles, cycles_edges = self.cycle_decomposition( fermion_edges_set )

        return (-1)**(len(cycles))

    @property
    def residue_type( self ):
        def dir_e(e, v):
            if self.edge_weights[e] == 2:   return 1
            if v == self.edges[e][0]:       return -1
            else:                           return 1

        ext_types = [ dir_e(e,v) * self.edge_weights[e] for v in self.external_vtcs_set for e in self.adj_edges( v, self.edges_set ) ]
        return tuple(sorted(ext_types))

    def get_vtx_type( self, v ):
        def dir1(e, v):
            if self.edge_weights[e] == 2:   return 1
            if v == self.edges[e][0]:       return -1
            else:                           return 1

        def dir2(e, v):
            if self.edge_weights[e] == 2:   return 1
            if v == self.edges[e][0]:       return 1
            else:                           return -1
        
        adj_types = [ dir1(e,v)*self.edge_weights[e] for e in self.adj_edges( v, self.edges_set ) ]
        adj_types += [ dir2(e,v)*self.edge_weights[e] for e in self.edges_set if self.edges[e] == (v,v) ]

        return tuple(sorted(adj_types))

    def get_vtcs_coloring( self ):
        dictWeights = { edge : self.edge_weights[e] for e,edge in enumerate(self.edges) }
        edge_degree_counter = self.edge_degree_counter(self.edges_set)
        selfloop_degree_list = [ (edge_degree_counter[(v,v)],dictWeights[(v,v)] if edge_degree_counter[(v,v)] else 2) for v in self.internal_vtcs_set ]

        selfloop_multiplicity_list = sorted( (mul,v) for v, mul in zip(self.internal_vtcs_set, selfloop_degree_list) )
        ( ( max_selfloop_multiplicity, _), _ ) = selfloop_multiplicity_list[-1] if selfloop_multiplicity_list else ((0,2), 0)

        self_loop_list = [ frozenset( vtx for mul, vtx in filter( lambda ((mul, we), vtx) : mul == i and we == w, selfloop_multiplicity_list ) ) for i in range( max_selfloop_multiplicity+1 ) for w in (1,2,3) ]
        
        return self_loop_list + [ self.external_vtcs_set ]

    def get_edges_coloring( self, edges_set ):
        fermion_edges_set = self.sub_edges_by_weight(1) & edges_set
        boson_edges_set = self.sub_edges_by_weight(2) & edges_set
        ghost_edges_set = self.sub_edges_by_weight(3) & edges_set

        fermion_edges = frozenset( self.edges[i] for i in fermion_edges_set if not self.is_selfloop(self.edges[i]) )
        ghost_edges = frozenset( self.edges[i] for i in ghost_edges_set if not self.is_selfloop(self.edges[i]) )
        boson_edges = frozenset( self.edges[i] for i in boson_edges_set )

        normalize = lambda edge : (max(edge),min(edge))
        flip = lambda (x,y) : (y,x)
        fermion_loops = frozenset( normalize(edge) for edge in fermion_edges if flip(edge) in fermion_edges )
        ghost_loops = frozenset( normalize(edge) for edge in ghost_edges if flip(edge) in ghost_edges )
        reduced_fermion_edges = fermion_edges - fermion_loops - frozenset( flip(edge) for edge in fermion_loops )
        reduced_ghost_edges = ghost_edges - ghost_loops - frozenset( flip(edge) for edge in ghost_loops )
        boson_fermion_loops = frozenset( edge for edge in reduced_fermion_edges if flip(edge) in boson_edges or edge in boson_edges )
        boson_ghost_loops = frozenset( edge for edge in reduced_ghost_edges if flip(edge) in boson_edges or edge in boson_edges )

        reduced_boson_edges = boson_edges - boson_fermion_loops - frozenset( flip(edge) for edge in boson_fermion_loops ) - boson_ghost_loops - frozenset( flip(edge) for edge in boson_ghost_loops )

        dbl_boson_edges = reduced_boson_edges | frozenset( flip(edge) for edge in reduced_boson_edges )

        if len(dbl_boson_edges&reduced_fermion_edges) != 0 or \
            len(dbl_boson_edges&reduced_ghost_edges) != 0:
            print dbl_boson_edges, reduced_fermion_edges
            raise

        boson_coloring = super( WeightedGraph, self).get_edges_coloring( boson_edges_set )

        return [ dbl_boson_edges | reduced_fermion_edges | reduced_ghost_edges, 
            fermion_loops, boson_fermion_loops, ghost_loops, boson_ghost_loops, 
            reduced_ghost_edges - boson_ghost_loops ] + boson_coloring[1:]

    def get_trivial_symmetry_factor( self ):
        grpSize = 1
        boson_edges = self.sub_edges_by_weight(2)
        edge_degree_counter = self.edge_degree_counter(boson_edges)
        for mul_edge_deg in ( m for edge, m in edge_degree_counter.iteritems() if not self.is_selfloop(edge) ):
            grpSize*= factorial(mul_edge_deg)

        for selfloop_deg in ( m for edge, m in edge_degree_counter.iteritems() if self.is_selfloop(edge) ):
            grpSize*= double_factorial(2*selfloop_deg)
        return grpSize

    def permute_external_edges( self ):
        class FixedGraph( type(self) ):
            def get_vtcs_coloring( self ):
                vtcs_coloring = super(FixedGraph, self).get_vtcs_coloring()

                vtcs_coloring = [ c - self.external_vtcs_set for c in vtcs_coloring]
                vtcs_coloring.extend( frozenset([v]) for v in sorted(self.external_vtcs_set) )

                return vtcs_coloring

        extern_boson_vtcs = \
            frozenset( v for e in self.sub_edges_by_weight(2) for v in self.edges[e] ) \
            & self.external_vtcs_set
        extern_in_fermion_vtcs = \
            frozenset( self.edges[e][0] for e in self.sub_edges_by_weight(1) ) \
            & self.external_vtcs_set
        extern_out_fermion_vtcs = \
            frozenset( self.edges[e][1] for e in self.sub_edges_by_weight(1) ) \
            & self.external_vtcs_set

        extern_vtcs_list =  list(extern_boson_vtcs) + \
                            list(extern_in_fermion_vtcs) + \
                            list(extern_out_fermion_vtcs)
        if frozenset(extern_vtcs_list) != self.external_vtcs_set:
            raise

        vtcs_list = list(self.internal_vtcs_set) + \
                    extern_vtcs_list

        for perm0 in itertools.permutations( extern_boson_vtcs ):
            for perm1 in itertools.permutations( extern_in_fermion_vtcs ):
                for perm2 in itertools.permutations( extern_out_fermion_vtcs ):

                    new_vtcs_list = tuple(self.internal_vtcs_set) + \
                                    perm0 + perm1 + perm2
                    m = dict( zip( vtcs_list, new_vtcs_list ) )

                    def relabel_edge( (v1,v2) ):
                        return (m[v1], m[v2])

                    yield FixedGraph( 
                            [ relabel_edge(edge) for edge in self.edges ], self.edge_weights, 0 )

    @property
    def clean_graph( self ):
        ext_sorter = ( e in self.external_edges_set for e,edge in enumerate(self.edges) )

        norm = lambda (edge) : (max(edge),min(edge))
        edges = [ norm(edge) if w == 2 else edge for w,edge in zip(self.edge_weights, self.edges) ]
        xwe_list = list(sorted(zip(ext_sorter, self.edge_weights, edges)))
        edges = [ edge for x,w,edge in xwe_list ]
        weights = [ w for x,w,edge in xwe_list ]
        g = copy.copy(self)
        g.edges = tuple(edges)
        g.edge_weights= tuple(weights)
        g.prepare_graph()

        return g
