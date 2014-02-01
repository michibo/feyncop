
import itertools
from graph import Graph
import collections

import nauty_ctrl

def calc_gen_params( num_loops, valence, num_ext_legs, cntd, notadpoles ):
    D = valence - 2

    num_vtcs = ( num_ext_legs + 2 * (num_loops -1 ) )
    num_edges = ( valence * ( num_loops -1 ) + num_ext_legs )

    if num_edges % D != 0 or num_vtcs % D != 0:
        return ( 0, 0, 0 )

    num_vtcs /= D
    num_edges /= D
    bulk_num_vtcs = num_vtcs + num_ext_legs
    bulk_num_edges = num_edges + num_ext_legs

    min_edges = 0
    max_edges = bulk_num_edges 

    if notadpoles:
        min_edges = bulk_num_edges
    elif cntd:
        min_edges = bulk_num_vtcs - 1
    elif valence % 2 == 1:
        min_edges = bulk_num_vtcs / 2
    else:
        min_edges = (num_ext_legs+1)/2

    return bulk_num_vtcs, min_edges, max_edges

def gen_graphs( num_loops, valence, num_ext_legs, cntd, edge2cntd, vtx2cntd, notadpoles ):
    cntd = cntd | edge2cntd | vtx2cntd
    edge2cntd = edge2cntd | vtx2cntd
    notadpoles = notadpoles | vtx2cntd

    bulk_num_vtcs, min_edges, max_edges = calc_gen_params( num_loops, valence, num_ext_legs, cntd, notadpoles )

    if bulk_num_vtcs <= 0:
        return

    for g_bulk in nauty_ctrl.gen_nauty_graphs( 
        bulk_num_vtcs, cntd, valence, min_edges, max_edges ):
        labeled_graphs = ( g for g in gen_from_bulk_g( 
                    g_bulk, frozenset(range(bulk_num_vtcs)), 
                    valence, num_ext_legs, notadpoles ) )

        unlabeled_graphs = frozenset( g.unlabeled_graph for g in labeled_graphs )

        for g in unlabeled_graphs:
            if vtx2cntd and not g.is_vtx_2_connected : continue
            elif edge2cntd and not g.is_edge_2_connected : continue
            elif cntd and not g.is_connected : continue

            if not vtx2cntd and notadpoles and g.is_tadpole : continue

            yield g

def gen_from_bulk_g( g, vtcs_set, valence, num_ext_legs, notadpoles ):
    valences = g.valency_dict
    leg_candidates = frozenset( v for v in vtcs_set if valences[v] == 1 )

    def gen_ext_vtcs():
        if len(leg_candidates) < num_ext_legs:
            return

        if len(leg_candidates) == num_ext_legs:
            yield leg_candidates
            return

        if not notadpoles and valence % 2 == 1: # Extra legs can still be closed with self-loops!
            for ext_vtcs in itertools.combinations(leg_candidates, num_ext_legs):
                yield frozenset(ext_vtcs)

    for ext_vtcs in gen_ext_vtcs():
        int_vtcs = vtcs_set - ext_vtcs

        degree_defs = [ valence - valences[v] for v in int_vtcs ]
        if any( d % 2 != 0 or d < 0 for d in degree_defs ):
            continue

        selfloop_edges = [ (v,v) for v,d in zip(int_vtcs, degree_defs) for i in range(d/2) ]

        edges = g.edges + selfloop_edges
        yield Graph( edges )

