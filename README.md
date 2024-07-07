Feyngen and Feyncop
===================


What is it? 
-----------

feyngen is a program to generate Feynman graphs for the use in perturbative 
calculations of quantum field theory.

feyncop is a program to calculate the coproduct of Feynman graphs in the 
scope of the Hopf algebra of Feynman graphs. 

Both programs are designed to be used at relatively high loop orders where 
traditional programs like [QGRAF](http://cfif.ist.utl.pt/~paulo/qgraf.html) are 
not applicable. 
For instance the established [nauty package](http://pallini.di.uniroma1.it/) is 
used to ensure high performance. 

The theoretical background to the programs with details to validation and implementation is outlined in my paper on [Feynman graph generation and calculations in the Hopf algebra of Feynman graphs][fcg] and in my master's thesis on [Algorithmization of the Hopf algebra of Feynman graphs](http://www2.mathematik.hu-berlin.de/~kreimer/wp-content/uploads/BorinskyMaster.pdf).

Please cite Michael Borinsky, [Feynman graph generation and calculations in the Hopf algebra of Feynman graphs][fcg], *Computer Physics Communications*, Volume 185, Issue 12, December 2014, Pages 3317–3330 if you want to refer to the programs. 
[fcg]: http://dx.doi.org/10.1016/j.cpc.2014.07.023

Download & Manual
-----------------

The source code for both programs can be downloaded from [github](https://github.com/michibo/feyncop). 
A [pre-built version](http://people.physik.hu-berlin.de/~borinsky/static/feyncop_built.tar.gz) and a seperate [manual](http://people.physik.hu-berlin.de/~borinsky/static/feyngencop_manual.pdf) for both programs are available on my webpage: http://people.physik.hu-berlin.de/~borinsky/

Prerequisites
-------------

To use either of the programs Python 3 with development files must be 
installed on your machine. For information on how to install Python please 
consult http://www.python.org/ .

Additionally if you do not use the pre-built version, the nauty package by Brendan McKay is needed. The newest version can be downloaded from: http://pallini.di.uniroma1.it/

Installation
------------

This step can be skipped if you want to use the pre-built version. 

Copy the feyncop and nauty archives into the same directory and extract them:

    $ tar xzf feyncop.tar.gz
    $ tar xzf nautyXXXX.tar.gz

Where XXXX are some letters representing the current nauty version. 

Next, change the name of the folder containing the nauty package:

    $ mv nautyXXXX/ nauty/

and build the nauty package:

    $ cd nauty/
    $ ./configure && make
    $ cd ../

Now, feyncop and feyngen can be build:

    $ cd feyncop/
    $ make

The two python programs feyngen and feyncop in the feyncop/ directory should 
know be working as expected. 

An overview of the parameters of the two programs is displayed with 

    $ ./feyngen --help

or

    $ ./feyncop --help

.

Testing
-------

To test feyngen run,

    $ ./feyngen 2 -j2

in the feyncop/ directory.

The output should be: 

    phi4_j2_h2 :=
    +G[[0,0],[0,0],[1,1],[1,1],[3,2]]/128
    +G[[0,0],[1,1],[1,1],[2,0],[3,0]]/16
    +G[[1,0],[1,0],[1,1],[2,0],[3,0]]/4
    +G[[0,0],[1,0],[1,0],[1,1],[3,2]]/16
    +G[[1,0],[1,0],[1,0],[1,0],[3,2]]/48
    +G[[0,0],[1,0],[1,1],[2,0],[3,1]]/4
    +G[[1,0],[1,0],[1,0],[2,0],[3,1]]/6
    ;

Corresponding to the sum of all 2-point, 2-loop diagrams in phi^4 
theory. 

To test feyncop run, 

    $ ./feyngen 2 -j2 -p | ./feyncop -u

the output should be:

    phi4_j2_h2_red_cop_unlab :=
    + 1/4 * T[ G[[0,0],[1,0],[2,0]], G[[0,0],[1,0],[2,0]] ]
    + 3/4 * T[ G[[1,0],[1,0],[2,0],[3,0],[4,1],[5,1]], G[[0,0],[1,0],[2,0]] ]
    ;

This output corresponds to the coproduct of the sum of all 1PI, 
2-point, 2-loop diagrams in phi^4 theory.


Included files
--------------

File                  | Description
----------------------|------------------------------------------------------
feyngen               | A program to generate Feynman graphs.
feyncop               | A program to calculate the coproduct of Feynman graphs.
graph.py              | Implements basic graph handling and algorithms.
weighted_graph.py     | Implements handling of QED and QCD graphs.
hopf_graph.py         | Implements the Hopf algebra properties of graphs.
phi_k_gen.py          | Code for the phi^k-theory graph generation.
phi_34_gen.py         | Code for the phi^3+phi^4-theory graph generation.
qed_gen.py            | Code for the QED graph generation.
qcd_gen.py            | Code for the QCD graph generation.
combinatorics.py      | Zero dimensional QFT calculations.
powerseries.py        | Ring of truncated power series calculations.
stuff.py              | Additional combinatorial helper functions.
nauty_ctrl.py         | Wrapper for the graph generation using geng and multig.
nauty_wrapper.c       | Wrapper code for the nauty canonical labeling function. 
Makefile              | Installation Makefile.
README                | This file.

