
"""qcd_gen.py: This file is part of the feyncop/feyngen package.
    Implements functions to generate QCD graphs. """

# See also: http://people.physik.hu-berlin.de/~borinsky/

__author__ = "Michael Borinsky"
__email__ = "borinsky@physik.hu-berlin.de"
__copyright__ = "Copyright (C) 2014 Michael Borinsky"
__license__ = "MIT License"
__version__ = "1.0"

# Copyright (c) 2014 Michael Borinsky

# This program is distributed under the MIT License:

# Permission is hereby granted, free of charge, to any person obtaining a copy 
# of this software and associated documentation files (the "Software"), to 
# deal in the Software without restriction, including without limitation the 
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or 
# sell copies of the Software, and to permit persons to whom the Software is 
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in 
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS 
# IN THE SOFTWARE.

import itertools
from weighted_graph import WeightedGraph

import phi_34_gen

def gen_graphs( L, r_t2, u_t2, m, cntd, edge2cntd, vtx2cntd, notadpoles ):
    """Generate QCD graphs with the desired parameters and properties.
        L: Loop number
        r_t2: Ext. fermion number
        u_t2: Ext. ghost number
        m: Ext. boson number"""
    
    phi34_graphs = ( phi_34_gen.gen_graphs( L, r_t2 + m + u_t2, cntd, edge2cntd, vtx2cntd, notadpoles) )

    for g_phi34 in phi34_graphs:
        gen_qcd_graphs = ( g.unlabeled_graph for g in gen_from_phi34_g( g_phi34, r_t2, u_t2, m ) )

        qcd_graphs = frozenset( gen_qcd_graphs )

        for g in qcd_graphs:
            yield g

def gen_from_phi34_g( fg, r_t2, u_t2, m ):
    """Helper function: Generate full fledged QCD graphs from the bulk output of 
        phi_34_gen.gen_graphs."""

    ext_vtcs = fg.external_vtcs_set
    int_vtcs = fg.internal_vtcs_set

    is_sl = [ fg.is_selfloop(edge) for e,edge in enumerate(fg.edges) ]
    ext_adj = [ frozenset(fg.adj_edges(v, fg.edges_set)) for v in ext_vtcs ]
    int_adj = [ frozenset(fg.adj_edges(v, fg.edges_set)) for v in int_vtcs ]
    for weights in itertools.product((1,2), repeat=len(fg.edges)):
        fermion_edges = frozenset( e for e,w in enumerate(weights) if w == 1 )
        fermion_adj = [ adj&fermion_edges for adj in int_adj ]
        fermion_valences = ( sum( 2 if is_sl[e] else 1 for e in adj ) for adj in fermion_adj )
        if any( (val != 0) and (val != 2) for val in fermion_valences ):
            continue

        boson_edges = frozenset( e for e,w in enumerate(weights) if w == 2 )
        boson_adj = [ adj&boson_edges for adj in int_adj ]
        boson_valences = ( sum( 2 if is_sl[e] else 1 for e in adj ) 
            for adj in boson_adj )
        if any( (val!=1) and (val!=3) and (val!=4) for val in boson_valences ):
            continue

        fermion_legs = sum( 1 for adj in ext_adj for e in adj if weights[e] == 1 )

        if fermion_legs != r_t2 + u_t2:
            continue

        boson_legs = sum( 1 for adj in ext_adj for e in adj if weights[e] == 2 )

        if boson_legs != m:
            continue

        for fermion_weights in itertools.product( (-1,1), repeat=len(fermion_edges)):
            dir_weights = list(weights)
            for i,e in enumerate(fermion_edges):
                dir_weights[e] = fermion_weights[i]
            
            def dir( e, v ):
                v1,v2 = fg.edges[e]
                return 1 if v1 == v else -1

            fermion_res = ( sum( 0 if is_sl[e] else dir(e,v)*dir_weights[e] for e in adj ) for v,adj in zip(int_vtcs,fermion_adj) )
            if any( res != 0 for res in fermion_res ):
                continue

            flip = lambda (v1,v2) : (v2,v1)
            edges = tuple( edge if w == 1 or w == 2 else flip(edge) for edge, w in zip(fg.edges, dir_weights) )
            translated_weights = tuple( 2 if w == 2 else 1 for w in weights )
            
            g = WeightedGraph( tuple(edges), tuple(translated_weights) )

            fermion_loops = list(g.cntd_components_sub_edges( g.sub_edges_by_weight(1) ))

            for ghost_loops in itertools.product( (True, False), repeat=len(fermion_loops) ):
                ghost_weights = list(g.edge_weights)
                for gh,loop in zip(ghost_loops, fermion_loops):
                    if gh:
                        for e in loop:
                            ghost_weights[e] = 3

                gw = WeightedGraph( tuple(edges), tuple(ghost_weights) )
                if len(gw.external_vtcs_set&gw.vtcs_set_sub_edges( gw.sub_edges_by_weight(3)&gw.external_edges_set)) == u_t2:
                    yield gw

