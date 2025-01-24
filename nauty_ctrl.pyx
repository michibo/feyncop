

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


def get_geng_obj(num_vtcs, cntd : str, max_degree, chunk=None):
    """
    Calls geng with the desired number of vertices, connectedness and
    maximum degree.

    The argument ``cntd`` can be ``"connected"`` or ``"biconnected"``
    or the empty string otherwise.

    The argument ``chunk`` should be either ``None`` or a pair of
    integers (i,N).  The integer N is a modulus. Only graphs with
    index equal to i mod N will be returned.
    """
    if cntd == "connected":
        cntd_param = ["-c"]
    elif cntd == "biconnected":
        cntd_param = ["-C"]
    else:
        cntd_param = []

    command = [str(nauty_path / "geng")]
    command.extend(cntd_param)
    command += ["-D%d" % max_degree, "-q", "%d" % num_vtcs]
    if chunk is not None:
        a, b = chunk
        a = int(a)
        b = int(b)
        command.append(f"{a}/{b}")
    return sp.Popen(command,
                    stdout=sp.PIPE, stderr=None, stdin=None)


def get_multig_obj(geng_stream, min_edges, max_edges, max_degree):
    """Calls multig with desired range of edges and max degree."""

    multig_path = str(nauty_path / "multig")
    return sp.Popen([multig_path, "-e%d:%d" % (min_edges, max_edges),
                     "-D%d" % max_degree, "-T", "-q"],
                    stdout=sp.PIPE, stderr=None, stdin=geng_stream)


cpdef multig_to_graph(multig_line):
    """
    Read the output of ``multig`` to a ``Graph`` object.

￼   EXAMPLES::

        sage: line = b'4 3  0 3 2 1 3 1 2 3 1\n'
￼       sage: multig_to_graph(line)
￼       G[[3,0],[3,0],[3,1],[3,2]]
    """
    cdef int v1, v2, mul

    cdef list info = multig_line.split()

    cdef list edge_info = info[2:]
    cdef list edges = []

    for x, y, z in zip(edge_info[::3], edge_info[1::3], edge_info[2::3]):
        v1, v2, mul = int(x), int(y), int(z)

        if v1 == v2:
            mul >>= 1

        edges.extend([(v2, v1)] * mul)

    return edges


def gen_nauty_graphs(num_vtcs, cntd : str, max_degree,
                     min_edges, max_edges, chunk=None):
    """
    Generate graphs with the desired properties.

    EXAMPLES::

        sage: from nauty_ctrl import *
        sage: list(gen_nauty_graphs(4,"connected",4,4,4))
        [G[[3,0],[3,0],[3,1],[3,2]],
         G[[2,0],[3,0],[3,0],[3,1]],
         G[[2,0],[2,0],[3,0],[3,1]],
         G[[2,0],[3,0],[3,1],[3,2]],
         G[[2,0],[3,0],[2,1],[3,1]]]
        sage: list(gen_nauty_graphs(4,"biconnected",4,4,4))
        [G[[2,0],[3,0],[2,1],[3,1]]]
    """
    geng = get_geng_obj(num_vtcs, cntd, max_degree, chunk=chunk)
    multig = get_multig_obj(geng.stdout, min_edges, max_edges, max_degree)

    for line in multig.stdout:
        yield Graph(multig_to_graph(line))


# SageMath alternative


def gen_nauty_graphs_sage(num_vtcs, cntd : str, max_degree,
                          min_edges, max_edges, chunk=None):
    """
    Generate SageMath graphs with the desired properties.

    EXAMPLES::

        sage: from nauty_ctrl import *
        sage: L = list(gen_nauty_graphs_sage(4,"connected",4,4,4)); L
        [Looped multi-graph on 4 vertices,
         Looped multi-graph on 4 vertices,
         Looped multi-graph on 4 vertices,
         Looped multi-graph on 4 vertices,
         Looped multi-graph on 4 vertices]
        sage: L[0].edges()
        [(0, 3, None), (0, 3, None), (1, 3, None), (2, 3, None)]
    """
    from sage.graphs.graph import Graph as SageGraph

    geng = get_geng_obj(num_vtcs, cntd, max_degree, chunk=chunk)
    multig = get_multig_obj(geng.stdout, min_edges, max_edges, max_degree)

    for line in multig.stdout:
        yield SageGraph(multig_to_graph(line), loops=True, multiedges=True,
                        format="list_of_edges")
