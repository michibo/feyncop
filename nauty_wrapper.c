
// nauty_wrapper.c: A wrapper for the nauty "nauty" routine to calculate 
// This file is part of the feyncop/feyngen package.
// the canonical labeling of a graph. 

// Author: Michael Borinsky
// Email: borinsky@physik.hu-berlin.de

#include <python3/Python.h>
#include "../nauty/nauty.h"

unsigned long long g_GroupSize;

static void grouplevelproc(int* p1,int* p2,int i1,int* p3,statsblk* s1,int i2,int index,int i3,int i4,int i5,int i6)
{
    g_GroupSize*= index;
}

static PyObject*
get_canonical_labeling(PyObject *self, PyObject *args)
{
    int num_vtcs = 0;
    PyObject *py_edges_list;
    PyObject *py_partition_list;

    if (!PyArg_ParseTuple(args, "iOO", &num_vtcs, &py_edges_list, &py_partition_list)) {
        PyErr_SetString(PyExc_TypeError, "Missing argument.");
        return NULL;
    }

    int n = num_vtcs;

    DYNALLSTAT(graph,g,g_sz);
    DYNALLSTAT(graph,cg,cg_sz);
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    DYNALLSTAT(setword,workspace,workspace_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    options.userlevelproc = grouplevelproc;


    options.getcanon = TRUE;
    options.defaultptn = FALSE; 
    options.digraph = TRUE;


    int m = (n + WORDSIZE - 1) / WORDSIZE;


    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);
    DYNALLOC2(graph, g, g_sz, m, n, "malloc");
    DYNALLOC2(graph, cg, cg_sz, m, n, "malloc");
    DYNALLOC1(setword, workspace, workspace_sz, 500*m, "malloc");
    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, n,"malloc");

    int* p_ptn = ptn;
    int* p_lab = lab;

    PyObject* partition_iter = PyObject_GetIter(py_partition_list);
    PyObject* partition;

    int num_vtcs2 = 0;
    while( (partition = PyIter_Next(partition_iter)) )
    {
        PyObject* vtx_iter = PyObject_GetIter(partition);

        PyObject* vtx;
        while( (vtx = PyIter_Next( vtx_iter )) )
        {
            int num = PyInt_AS_LONG( vtx );
            Py_DECREF( vtx );

            *(p_lab++) = num;
            *(p_ptn++) = 1;

            num_vtcs2++;
        }

        *(p_ptn-1) = 0;

        Py_DECREF( vtx_iter );
        Py_DECREF( partition );
    }
    Py_DECREF( partition_iter );

    if ( num_vtcs != num_vtcs2 )
    {
        PyErr_SetString(PyExc_TypeError,
                "vtx number doesn't fit to partition");
        return 0;       // error
    }

    int v;
    for ( v = 0; v < n; v++ )
    {
        set *gv = GRAPHROW(g,v,m);
        EMPTYSET(gv,m);
    }

    PyObject* edges_iter = PyObject_GetIter( py_edges_list );

    PyObject* edge;
    while( (edge = PyIter_Next(edges_iter)) )
    {
        PyObject* v1 = PyTuple_GetItem( edge, 0 );
        PyObject* v2 = PyTuple_GetItem( edge, 1 );
        int iv1 = PyInt_AS_LONG( v1 );
        int iv2 = PyInt_AS_LONG( v2 );
        Py_DECREF( edge );

        ADDELEMENT((GRAPHROW(g, iv1, m)), iv2);
    }
    Py_DECREF( edges_iter );

    g_GroupSize = 1;

    // Calling nauty.
    nauty(g, lab, ptn, NULL, orbits, &options, &stats, workspace, 500*m, m, n, cg);
 
    int i;

    // Reading the canonical labeling and the orbits.
    PyObject* py_lab= PyList_New(num_vtcs);
    for ( i=0; i < num_vtcs; i++ )
    {
        PyList_SET_ITEM( py_lab, i, Py_BuildValue("i", (int)lab[i]) );
    }

    PyObject* py_orbits = PyList_New(num_vtcs);
    for( i=0; i < num_vtcs; i++ )
    {
        PyList_SET_ITEM( py_orbits, i, Py_BuildValue("i", (int)orbits[i]) );
    }
    DYNFREE(g, g_sz);
    DYNFREE(cg, cg_sz);
    DYNFREE(workspace, workspace_sz);
    DYNFREE(lab, lab_sz);
    DYNFREE(ptn, ptn_sz);
    DYNFREE(orbits, orbits_sz);

    PyObject* pyret = PyTuple_New(3);
    PyTuple_SET_ITEM( pyret, 0, py_lab );
    PyTuple_SET_ITEM( pyret, 1, Py_BuildValue("K", g_GroupSize ) );
    PyTuple_SET_ITEM( pyret, 2, py_orbits );

    return pyret;
}

static PyMethodDef nauty_wrapper_funcs[] = {
    {"get_canonical_labeling", get_canonical_labeling, METH_VARARGS, "Give canonical labeling" },
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC
initnauty_wrapper(void) {
    Py_InitModule( "nauty_wrapper", nauty_wrapper_funcs );
}

