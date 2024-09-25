#!/usr/bin/env python

""" - Main program code for "feyngen" - A program to generate Feynman graphs
suitable for Hopf algebra calculations. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

import argparse
import sys

from sage.all import factorial, QQ, CombinatorialFreeModule, NN, Words
import combinatorics
from weighted_graph import WeightedGraph
import phi_k_gen
import qed_gen
import qcd_gen
import phi_34_gen


def gen_and_count(gen_graphs, non_leg_fixed, test_sums):
    """Generates the graphs and sums the symmetry factors for the obligatory
        cross check with the zero-dimensional field theory terms."""

    sym_fac_non_fixed = QQ((0, 1))
    sym_fac_fixed = QQ((0, 1))

    for g in gen_graphs:
        sym_fac_non_fixed += QQ((1, g.symmetry_factor))

        if not non_leg_fixed:
            gen_fixed_graphs = (fixed_g.unlabeled_graph for fixed_g in g.permute_external_edges())

            fixed_graphs = frozenset(gen_fixed_graphs)
            for g_fixed in fixed_graphs:
                sym_fac_fixed += QQ((1, g_fixed.symmetry_factor))
                yield g_fixed
        else:
            yield g

    test_sums['fixed'] = sym_fac_fixed
    test_sums['non-fixed'] = sym_fac_non_fixed


def generating_phi_k(loops, ext_legs, valence=4, connected=False):
    """
    sketch of generator for the class Phi_k
    """
    gen_graphs = phi_k_gen.gen_graphs(loops,
                                      valence,
                                      ext_legs,
                                      connected,
                                      False,
                                      False,
                                      False)
    test_sums = dict()
    return gen_and_count(gen_graphs, False, test_sums)


def sum_over_graphs(graphs):
    """
    Return the sum over graphs.

    EXAMPLES::

        sage: from feyngen_sage import *
        sage: L = generating_phi_k(2, 2)
        sage: sum_over_graphs(L)
        1/128*B[word: (0, 0),(0, 0),(1, 1),(1, 1),(2, 3)] + 1/16*B[word: (0, 0),(0, 1),(0, 1),(1, 1),(2, 3)] + 1/4*B[word: (0, 0),(0, 1),(0, 2),(1, 1),(1, 3)] + 1/16*B[word: (0, 0),(0, 2),(0, 3),(1, 1),(1, 1)] + 1/48*B[word: (0, 1),(0, 1),(0, 1),(0, 1),(2, 3)] + 1/6*B[word: (0, 1),(0, 1),(0, 1),(0, 2),(1, 3)] + 1/4*B[word: (0, 1),(0, 1),(0, 2),(0, 3),(1, 1)]
        sage: sum(_.coefficients())
    """
    indices = Words(NN.cartesian_product(NN), infinite=False)
    M = CombinatorialFreeModule(QQ, indices)

    def short_str(monome):
        return "G" + str(tuple(monome))

    M._repr_term = short_str
    word = M.basis().keys()

    return M.sum_of_terms((word(g.get_edges_tuple()),
                           QQ.one() / g.symmetry_factor)
                          for g in graphs)


def compare_sym_factors(num_loops, test_sums, args):
    """Compares the sum of the symmetry factors to the corresponding
        combinatorial calculations using generating functions.
        Returns False if a discrepancy is found."""

    test_sum_a = 0
    fix_factor = 0

    if args.ym:
        fix_factor = factorial(args.num_ext_blegs) * factorial(args.num_ext_flegs // 2)**2 * factorial(args.num_ext_glegs // 2)**2
        if args.connected:
            test_sum_a = combinatorics.cntd_qcd_class_coeff(num_loops, args.num_ext_flegs, args.num_ext_glegs, args.num_ext_blegs)
        else:
            test_sum_a = combinatorics.qcd_class_coeff(num_loops, args.num_ext_flegs, args.num_ext_glegs, args.num_ext_blegs)
    elif args.qed_furry:
        fix_factor = factorial(args.num_ext_blegs) * factorial(args.num_ext_flegs // 2)**2
        if args.connected:
            test_sum_a = combinatorics.cntd_qed_furry_class_coeff(num_loops, args.num_ext_flegs, args.num_ext_blegs)
        else:
            test_sum_a = combinatorics.qed_furry_class_coeff(num_loops, args.num_ext_flegs, args.num_ext_blegs)
    elif args.qed:
        fix_factor = factorial(args.num_ext_blegs) * factorial(args.num_ext_flegs // 2)**2
        if args.connected:
            test_sum_a = combinatorics.cntd_qed_class_coeff(num_loops, args.num_ext_flegs, args.num_ext_blegs)
        else:
            test_sum_a = combinatorics.qed_class_coeff(num_loops, args.num_ext_flegs, args.num_ext_blegs)
    elif args.phi34:
        fix_factor = factorial(args.num_ext_legs)
        if args.connected:
            test_sum_a = combinatorics.cntd_phi34_class_coeff(num_loops, args.num_ext_legs)
        else:
            test_sum_a = combinatorics.phi34_class_coeff(num_loops, args.num_ext_legs)
    else:
        fix_factor = factorial(args.num_ext_legs)
        if args.connected:
            test_sum_a = combinatorics.cntd_phi_k_class_coeff(num_loops, args.num_ext_legs, args.valence)
        else:
            test_sum_a = combinatorics.phi_k_class_coeff(num_loops, args.num_ext_legs, args.valence)

    if test_sums['non-fixed'] != test_sum_a:
        print(test_sums['non-fixed'], test_sum_a)
        print('nonfixed')
        return False

    if not args.non_leg_fixed and test_sums['fixed'] != test_sum_a * fix_factor:
        print(test_sums['fixed'], test_sum_a * fix_factor)
        print('fixed', fix_factor)
        return False

    return True


def main():
    """Main program section. Reads the options and parameters and starts the
        relevant subroutines."""

    parser = argparse.ArgumentParser(description='Generate non-isomorphic feynman diagrams with their corresponding symmetry-factors.')
    parser.add_argument('loops', metavar='#loops', type=int, nargs='+', help='Loop number(s) to generate or order in h-bar for non connected graphs')

    parser.add_argument('-c', '--connected', dest='connected', action='store_true', help='Generate only connected graphs. (default: false)')
    parser.add_argument('-p', '--1PI', dest='edge2cntd', action='store_true', help='Generate only 1PI graphs. Implies -c/--connected. (default: false)')
    parser.add_argument('-v', '--vtx2cntd', dest='vtx2cntd', action='store_true', help='Generate only 2-vertex connected graphs. Implies -p/--1PI and -t/--notadpoles. (default: false)')
    parser.add_argument('-t', '--notadpoles', dest='notadpoles', action='store_true', help='Generate only non-tadpole graphs. (default: false)')

    parser.add_argument('-k', '--valence', dest='valence', type=int, default=4, help='Generate phi^k graphs. k determines the vertex valence of the phi^k-QFT. (default: 4)')
    parser.add_argument('--phi34', dest='phi34', action='store_true', help='Generate graphs with only 3 or 4 valent internal vertices. -k is ignored. (default: false)')

    parser.add_argument('--qed', dest='qed', action='store_true', help='Generate QED graphs (with neglect of Furry\'s theorem). -k is ignored. (default: false)')
    parser.add_argument('--qed_furry', dest='qed_furry', action='store_true', help='Generate QED graphs (respecting Furry\'s theorem). -k is ignored. (default: false)')
    parser.add_argument('--ym', dest='ym', action='store_true', help='Generate Yang-Mills graphs. -k is ignored. (default: false)')

    parser.add_argument('-j', '--ext_legs', dest='num_ext_legs', type=int, default=0, help='Set the number of external legs in scalar phi^k-QFT. (default: 0)')
    parser.add_argument('-b', '--ext_boson_legs', dest='num_ext_blegs', type=int, default=0, help='Set the number of external photon/gluon legs (default: 0)')
    parser.add_argument('-f', '--ext_fermion_legs', dest='num_ext_flegs', type=int, default=0, help='Set the number of external fermion legs (default: 0)')
    parser.add_argument('-g', '--ext_ghost_legs', dest='num_ext_glegs', type=int, default=0, help='Set the number of external ghost legs (default: 0)')

    parser.add_argument('-u', '--non_leg_fixed', dest='non_leg_fixed', action='store_true', help='Generate non-leg-fixed graphs. I.e. don\'t label/distinguish external legs of graphs during isomorphism testing and symmetry-factor calculation. (default: false)')

    args = parser.parse_args()

    qed = args.qed or args.qed_furry

    if qed or args.ym:
        if args.num_ext_legs:
            print("Warning: -j# is ignored. Generating QED graphs, use -b and -f for # of boson and fermion legs instead.")
    if args.num_ext_legs and (args.num_ext_blegs or args.num_ext_flegs or args.num_ext_glegs):
        print("Warning: -b# and -f# are ignored. Use only -b, -f or -g with --qed or --ym and without -j to indicate QED or Yang Mills enumeration.")
    elif (args.num_ext_blegs or args.num_ext_flegs) and not args.ym:
        args.qed = qed = True

    if not qed and args.valence <= 2:
        print("Sorry, this program is designed to generate feynman graphs. The vertex valence (-k) needs to be >=3")
        return

    loops_str = ""
    loops_str = "%s" % "_".join(["h%d" % l for l in args.loops])

    unlabeled_str = "_nlf" if args.non_leg_fixed else ""
    if args.qed_furry:
        print("qed_with_furry_f%d_b%d%s_%s :=" % (args.num_ext_flegs, args.num_ext_blegs, unlabeled_str, loops_str))
    elif qed:
        print("qed_f%d_b%d%s_%s :=" % (args.num_ext_flegs, args.num_ext_blegs, unlabeled_str, loops_str))
    elif args.ym:
        print("ym_f%d_g%d_b%d%s_%s :=" % (args.num_ext_flegs, args.num_ext_glegs, args.num_ext_blegs, unlabeled_str, loops_str))
    elif args.phi34:
        print("phi3_4_j%d%s_%s :=" % (args.num_ext_legs, unlabeled_str, loops_str))
    else:
        print("phi%d_j%d%s_%s :=" % (args.valence, args.num_ext_legs, unlabeled_str, loops_str))

    ext_legs_total = 0
    possible_emptyness = not args.connected and not args.edge2cntd and not args.vtx2cntd

    errors_encountered = False

    for num_loops in args.loops:
        gen_graphs = None

        if qed:
            ext_legs_total = args.num_ext_flegs + args.num_ext_blegs
            gen_graphs = qed_gen.gen_graphs(num_loops, args.num_ext_flegs, args.num_ext_blegs, args.connected, args.edge2cntd, args.vtx2cntd, args.notadpoles, args.qed_furry)
        elif args.ym:
            ext_legs_total = args.num_ext_flegs + args.num_ext_blegs + args.num_ext_glegs
            gen_graphs = qcd_gen.gen_graphs(num_loops, args.num_ext_flegs, args.num_ext_glegs, args.num_ext_blegs, args.connected, args.edge2cntd, args.vtx2cntd, args.notadpoles)
        elif args.phi34:
            ext_legs_total = args.num_ext_legs
            gen_graphs = phi_34_gen.gen_graphs(num_loops, args.num_ext_legs, args.connected, args.edge2cntd, args.vtx2cntd, args.notadpoles)
        else:
            ext_legs_total = args.num_ext_legs
            gen_graphs = phi_k_gen.gen_graphs(num_loops, args.valence, args.num_ext_legs, args.connected, args.edge2cntd, args.vtx2cntd, args.notadpoles)

        test_sums = dict()
        gen_graphs_c = gen_and_count(gen_graphs, args.non_leg_fixed, test_sums)

        if sum_over_graphs(gen_graphs_c):

            testable = not args.notadpoles and not args.edge2cntd and not args.vtx2cntd
            if testable and not compare_sym_factors(num_loops, test_sums, args):
                errors_encountered |= True
                print("Warning: Internal error check failed", file=sys.stderr)

        elif num_loops == 1 and ext_legs_total == 0 and possible_emptyness:
            print("+%s" % str(WeightedGraph(tuple(), tuple()).unlabeled_graph))

            test_sums['fixed'] = test_sums['non-fixed'] = 1
            if not compare_sym_factors(num_loops, test_sums, args):
                errors_encountered |= True
                print("Warning: Internal error check failed", file=sys.stderr)
