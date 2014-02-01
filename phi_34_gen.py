
import itertools
from graph import Graph
import collections

import nauty_ctrl

def calc_gen_params( num_loops, num_ext_legs, cntd, notadpoles ):
    min_num_vtcs = ( num_ext_legs + 2*(num_loops - 1) + 1 ) / 2 + num_ext_legs
    max_num_vtcs = ( num_ext_legs + 2*(num_loops - 1) ) + num_ext_legs

    min_num_edges = ( 4*(num_loops - 1) + num_ext_legs + 1) / 2 + num_ext_legs/2
    max_num_edges = 3*(num_loops - 1) + 2*num_ext_legs

    if not notadpoles:
        if cntd:
            min_num_edges = min_num_vtcs - 1
        else:
            min_num_edges = (num_ext_legs+1)/2

    return min_num_vtcs, max_num_vtcs, min_num_edges, max_num_edges

def gen_graphs( num_loops, num_ext_legs, cntd, edge2cntd, vtx2cntd, notadpoles ):
    cntd = cntd | edge2cntd | vtx2cntd
    edge2cntd = edge2cntd | vtx2cntd
    notadpoles = notadpoles | vtx2cntd

    min_num_vtcs, max_num_vtcs, min_edges, max_edges = calc_gen_params( num_loops, num_ext_legs, cntd, notadpoles )

    if max_num_vtcs <= 0:
        return

    for bulk_num_vtcs in range(min_num_vtcs, max_num_vtcs+1):
        for g_bulk in nauty_ctrl.gen_nauty_graphs( 
            bulk_num_vtcs, cntd, 4, min_edges, max_edges ):
            labeled_graphs = ( g for g in gen_from_bulk_g( 
                        g_bulk, frozenset(range(bulk_num_vtcs)), 
                        num_loops, num_ext_legs, notadpoles ) )

            unlabeled_graphs = frozenset( g.unlabeled_graph for g in labeled_graphs )

            for g in unlabeled_graphs:
                if vtx2cntd and not g.is_vtx_2_connected : continue
                elif edge2cntd and not g.is_edge_2_connected : continue
                elif cntd and not g.is_connected : continue

                if not vtx2cntd and notadpoles and g.is_tadpole : continue

                yield g

def gen_from_bulk_g( g, vtcs_set, num_loops, num_ext_legs, notadpoles ):
    valences = g.valency_dict
    leg_candidates = frozenset( v for v in vtcs_set if valences[v] == 1 )

    def gen_ext_vtcs():
        if len(leg_candidates) < num_ext_legs:
            return

        if len(leg_candidates) == num_ext_legs:
            yield leg_candidates
            return

        if not notadpoles: # Extra legs can still be closed with self-loops!
            for ext_vtcs in itertools.combinations(leg_candidates, num_ext_legs):
                yield frozenset(ext_vtcs)

    for ext_vtcs in gen_ext_vtcs():
        int_vtcs = vtcs_set - ext_vtcs

        degree_defs = [ 4 - valences[v] for v in int_vtcs ]
        if any( d < 0 for d in degree_defs ):
            continue

        selfloop_edges = [ (v,v) for v,d in zip(int_vtcs, degree_defs) for i in range(d/2) ]

        edges = g.edges + selfloop_edges
        if len( edges ) - len(vtcs_set) == num_loops - 1:
            yield Graph( edges )

