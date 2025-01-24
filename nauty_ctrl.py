

"""nauty_ctrl.py: This file is part of the feyncop/feyngen package.
    Wrapper for the nauty programs geng and multig. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

import subprocess as sp
from pathlib import Path

from graph import Graph


nauty_path = Path(__file__).resolve().parent


def get_geng_obj(num_vtcs, cntd, max_degree):
    """Calls geng with the desired number of vertices, connectedness and
    maximum degree."""

    cntd_param = ["-c"] if cntd else []

    geng_path = str(nauty_path / "geng")
    return sp.Popen([geng_path, "%d" % num_vtcs,
                     "-D%d" % max_degree, "-q"] + cntd_param,
                    stdout=sp.PIPE, stderr=None, stdin=None)


def get_multig_obj(geng_stream, min_edges, max_edges, max_degree):
    """Calls multig with desired range of edges and max degree."""

    multig_path = str(nauty_path / "multig")
    return sp.Popen([multig_path, "-e%d:%d" % (min_edges, max_edges),
                     "-D%d" % max_degree, "-T", "-q"],
                    stdout=sp.PIPE, stderr=None, stdin=geng_stream)


def multig_to_graph(multig_line):
    """Reads the output of multig to a Graph object."""

    info = multig_line.split()

    num_edges = int(info[1])

    edge_info = info[2:3 * num_edges + 2]
    edges = []

    for x, y, z in zip(edge_info[::3], edge_info[1::3], edge_info[2::3]):
        v1, v2, mul = int(x), int(y), int(z)

        if v1 == v2:
            mul >>= 1

        edges.extend([(v2, v1)] * mul)

    return edges


def gen_nauty_graphs(num_vtcs, cntd, max_degree, min_edges, max_edges):
    """
    Generate graphs with the desired properties.

    EXAMPLES::

        sage: from nauty_ctrl import *
        sage: list(gen_nauty_graphs(4,True,4,4,4))
        [G[[3,0],[3,0],[3,1],[3,2]],
         G[[2,0],[3,0],[3,0],[3,1]],
         G[[2,0],[2,0],[3,0],[3,1]],
         G[[2,0],[3,0],[3,1],[3,2]],
         G[[2,0],[3,0],[2,1],[3,1]]]
    """
    geng = get_geng_obj(num_vtcs, cntd, max_degree)
    multig = get_multig_obj(geng.stdout, min_edges, max_edges, max_degree)

    for line in multig.stdout:
        yield Graph(multig_to_graph(line))


# SageMath alternative


def gen_nauty_graphs_sage(num_vtcs, cntd, max_degree, min_edges, max_edges):
    """
    Generate SageMath graphs with the desired properties.

    EXAMPLES::

        sage: from nauty_ctrl import *
        sage: L = list(gen_nauty_graphs_sage(4,True,4,4,4)); L
        [Looped multi-graph on 4 vertices,
         Looped multi-graph on 4 vertices,
         Looped multi-graph on 4 vertices,
         Looped multi-graph on 4 vertices,
         Looped multi-graph on 4 vertices]
        sage: L[0].edges()
        [(0, 3, None), (0, 3, None), (1, 3, None), (2, 3, None)]
    """
    from sage.graphs.graph import Graph as SageGraph

    geng = get_geng_obj(num_vtcs, cntd, max_degree)
    multig = get_multig_obj(geng.stdout, min_edges, max_edges, max_degree)

    for line in multig.stdout:
        yield SageGraph(multig_to_graph(line), loops=True, multiedges=True,
                        format="list_of_edges")
