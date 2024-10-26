
""" parsefg.py: This file is part of the feyncop/feyngen package.
    Implements the Graph parsing. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

import re
import sys
from fractions import Fraction
import collections

from hopf_graph import HopfGraph

var_name_ptrn = re.compile(r"^(\S+)\s*:=\s*")


def parse_var_name(string):
    """Parses the name of a input graph sum."""

    m = var_name_ptrn.search(string)

    if m:
        return m.group(1), m.end(1), m.end()
    else:
        return "", 0, 0


fraction_pattern = re.compile(r"^(\d*)(/?)(\d*)$")


def parse_fraction(s):
    """Parses a fraction encoded in a string s."""

    if not s:
        return Fraction(1, 1)

    m = fraction_pattern.match(s)

    n = 1
    try:
        n = int(m.group(1))
    except ValueError:
        pass

    d = 1
    try:
        d = int(m.group(3))
    except ValueError:
        pass

    return Fraction(n, d)


edge_pattern = re.compile(r"\[\s*(\d+)\s*,\s*(\d+)\s*(,?)\s*(-?\s*[Afc]+|)\s*\]")
graph_pattern = re.compile(r"\s*(\+?-?)\s*(\d*/?\d*)\s*\*?\s*G\[([0-9,\[\]\sAfc]*)\]\*?(\d*/?\d*)?\s*")
graph_pattern_pow = re.compile(r"\s*\(?(G\[[0-9,\[\]\sAfc]*\])\s*\)?\^?(\d*)")
tensor_product_pattern = re.compile(r"\s*(\+?-?)\s*(\d*/?\d*)\s*\*?\s*T\[\s*((?:\(?G\[[0-9,\[\]\sAfc]*\]\)?\^?(\d*)\s*\*?\s*)+),\s*(G\[[0-9,\[\]\sAfc]*\])\s*\]\s*")
graph_with_tp_pattern = re.compile(r"\s*(\+?-?)\s*(\d*/?\d*)\s*\*?\s*(G\[[0-9,\[\]\sAfc]*\])\s*\*\s*\(((?:\s*\+?-?\s*\d*/?\d*\s*\*?\s*T\[\s*(?:\(?(?:G\[[0-9,\[\]\sAfc]*\])\)?\^?\d*\s*\*?\s*)+\s*,\s*G\[[0-9,\[\]\sAfc]*\]\s*\])*)\s*\)\s*")


def get_graph_from_match(m):
    """Helper function: Parses a graph from a match."""

    edges_string = m.group(3)

    dict_W = {'A': 2, 'f': 1, 'c': 3}
    global ym
    ym = False

    def gen_edges():
        for m_e in edge_pattern.finditer(edges_string):
            v1 = int(m_e.group(1))
            v2 = int(m_e.group(2))
            w = dict_W[m_e.group(4)] if m_e.group(3) == "," else 2

            if m_e.group(3) == ",":
                global ym
                ym = True

            yield (v1, v2, w)

    edges_weights = tuple(gen_edges())

    edges = [(v1, v2) for v1, v2, w in edges_weights]
    weights = [w for v1, v2, w in edges_weights]

    f1 = parse_fraction(m.group(2))
    f2 = parse_fraction(m.group(4))
    sign = -1 if "-" in m.group(1) else 1

    return HopfGraph(edges, weights, 0), sign * f1 * f2, ym


def get_tensor_product_from_match(m):
    """Helper function: Parses a tensor product from a match."""

    gprs = m.groups()

    f1 = parse_fraction(gprs[1])
    res_str = gprs[-1]
    res_graph, res_fac, res_ym = get_graph_from_match(graph_pattern.match(res_str))
    if res_fac != 1:
        print(f"Warning strange input: {m.group(0)}", file=sys.stderr)
        return None

    def gen_sgs():
        sbgrs_str = gprs[2]
        for sg_m in graph_pattern_pow.finditer(sbgrs_str):
            sg_str = sg_m.group(1)
            exp_str = sg_m.group(2)
            sg, sg_fac, sg_ym = get_graph_from_match(graph_pattern.match(sg_str))
            if sg_fac != 1 or sg_ym != res_ym:
                print(f"Warning strange input: {m.group(0)}", file=sys.stderr)
                continue

            p = 1 if not exp_str else int(exp_str)

            yield sg, p

    sgs = collections.Counter(dict(gen_sgs()))
    return (tuple(sorted(sgs.items())), res_graph), f1, ym


def get_graph_with_tp_from_match(m):
    gprs = m.groups()
    f1 = parse_fraction(gprs[1])

    g, g_fac, g_ym = get_graph_from_match(graph_pattern.match(gprs[2]))
    if g_fac != 1:
        print(f"Warning strange input: {m.group(0)}", file=sys.stderr)

    tps_str = gprs[3]

    def gen_tps():
        for tp_m in tensor_product_pattern.finditer(tps_str):
            if not tp_m:
                print(f"Warning strange input: {m.group(0)}", file=sys.stderr)
                continue
            tp, fac, ym = get_tensor_product_from_match(tp_m)

            if ym != g_ym:
                print(f"Warning strange input: {m.group(0)}", file=sys.stderr)
                continue

            yield tp, fac

    tp_sum = collections.Counter(dict(gen_tps()) if tps_str != (None,) else {})
    return g, tp_sum, f1, g_ym


def parse_sum_of_graphs(string):
    """Parses a graph sum."""
    for m in graph_pattern.finditer(string):
        yield get_graph_from_match(m), m.start(), m.end()


def parse_sum_of_tensor_products(string):
    """Parses a tensor product sum."""
    for m in tensor_product_pattern.finditer(string):
        yield get_tensor_product_from_match(m), m.start(), m.end()


def parse_sum_of_graph_with_tp(string):
    """Parses a graph with tensor product sum."""
    for m in graph_with_tp_pattern.finditer(string):
        yield get_graph_with_tp_from_match(m), m.start(), m.end()


def not_parsable_check(s):
    if s:
        print("\n********************************", file=sys.stderr)
        print(f"Warning: Could not parse this: {s}", file=sys.stderr)
        print("********************************", file=sys.stderr)


end_pattern = re.compile(r"0?;\s*(.*)$")


def parse_input_lines(instream, outstream, string, parser_fun=parse_sum_of_graphs):
    """Parses a stream of input."""

    for line in instream:
        string += line
        oldend = 0
        for g_fac, strbeg, strend in parser_fun(string):
            not_parsable_check(string[oldend:strbeg])
            oldend = strend

            yield g_fac
        string = string[oldend:]

    oldend = 0
    for g_fac, strbeg, strend in parser_fun(string):
        not_parsable_check(string[oldend:strbeg])
        oldend = strend

        yield g_fac
    string = string[oldend:]

    if string:
        m = end_pattern.match(string)
        if not m:
            not_parsable_check(string)
        else:
            not_parsable_check(m.group(1))
