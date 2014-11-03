

"""nauty_ctrl.py: This file is part of the feyncop/feyngen package.
    Wrapper for the nauty programs geng and multig. """

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

import subprocess as sp
import os
from graph import Graph

def get_geng_obj( num_vtcs, cntd, max_degree ):
    """Calls geng with the desired number of vertices, connectedness and maximum degree."""

    cntd_param = ["-c"] if cntd else []

    geng_path = os.path.join( os.path.dirname(os.path.realpath(__file__)), "geng" )
    geng_obj = sp.Popen([ geng_path, "%d" % num_vtcs, "-D%d" % max_degree, "-q" ] + cntd_param, stdout=sp.PIPE, stderr=None, stdin=None)
    return geng_obj

def get_multig_obj( geng_stream, min_edges, max_edges, max_degree ):
    """Calls multig with desired range of edges and max degree."""

    multig_path = os.path.join( os.path.dirname(os.path.realpath(__file__)), "multig" )
    multig_obj = sp.Popen([ multig_path, "-e%d:%d" % (min_edges, max_edges), "-D%d" % max_degree, "-T", "-q" ], stdout=sp.PIPE, stderr=None, stdin=geng_stream)
    return multig_obj

def multig_to_graph( multig_line ):
    """Reads the output of multig to a Graph object."""

    info = multig_line.split()

    num_vtcs = int(info[0])
    num_edges = int(info[1])

    edge_info = info[2:3*num_edges+2]
    edges = []

    for x,y,z in zip( edge_info[::3], edge_info[1::3], edge_info[2::3] ):
        v1,v2,mul = int(x), int(y), int(z)
        
        if v1 == v2:
            mul /= 2

        edges.extend( [(v2,v1)] * mul )

    return Graph(edges)

def gen_nauty_graphs( num_vtcs, cntd, max_degree, min_edges, max_edges ):
    """Generate graphs with the desired properties."""

    geng = get_geng_obj( num_vtcs, cntd, max_degree )
    multig = get_multig_obj( geng.stdout, min_edges, max_edges, max_degree )

    for line in multig.stdout:
        yield multig_to_graph( line )
