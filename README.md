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

The theoretical background to the programs with details to validation and implementation is outlined in the paper on [Feynman graph generation and calculations in the Hopf algebra of Feynman graphs](http://dx.doi.org/10.1016/j.cpc.2014.07.023).

Please cite Michael Borinsky, [Feynman graph generation and calculations in the Hopf algebra of Feynman graphs](http://dx.doi.org/10.1016/j.cpc.2014.07.023), *Computer Physics Communications*, Volume 185, Issue 12, December 2014, Pages 3317–3330 if you want to refer to the programs.

Acknowledgements
----------------

Many thanks are owed to [Frédéric Chapoton](//irma.math.unistra.fr/~chapoton/) for the long overdue upgrade of feyngen & feyncop to Python3 in 2024.

Bug reports or (pull) requests are always welcome.

Download & Manual
-----------------

The source code for both programs can be downloaded from [github](https://github.com/michibo/feyncop). 
A [deprecated python2 pre-built version](https://michaelborinsky.com/static/feyncop_built.tar.gz) and a separate (also slightly outdated) [manual](https://michaelborinsky.com/static/feyngencop_manual.pdf) for both programs are available on my webpage: https://michaelborinsky.com

Prerequisites
-------------

To use either of the programs Python 3 with development files must be
installed on your machine. For information on how to install Python please
consult https://www.python.org/

You need also a C compiler and the Cython compiler, so that both commands `gcc` and `cythonize` are available on your machine.

Additionally if you do not use the pre-built version, the nauty package by Brendan McKay is needed. The newest version can be downloaded from: http://pallini.di.uniroma1.it/

Installation
------------

This step can be skipped if you want to use the pre-built version.

Clone the feyncop repository and copy the nauty archives into the same directory and extract them:

    $ git clone https://github.com/michibo/feyncop.git
    $ tar xzf nautyXXXX.tar.gz

Where XXXX are some letters representing the current `nauty` version that can be downloaded from http://pallini.di.uniroma1.it/

Next, change the name of the folder containing the `nauty` package:

    $ mv nautyXXXX/ nauty/

and build the nauty package:

    $ cd nauty/
    $ ./configure && make
    $ cd ../

Now, feyncop and feyngen can be build:

    $ cd feyncop/
    $ make

The two python programs `feyngen` and `feyncop` in the `feyncop/` directory should
now be working as expected.

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

To quickly test feyncop run,

    $ ./feyngen 2 -j2 -p | ./feyncop -u

the output should be:

    phi4_j2_h2_red_cop_unlab :=
    + 1/4 * T[ G[[0,0],[1,0],[2,0]], G[[0,0],[1,0],[2,0]] ]
    + 3/4 * T[ G[[1,0],[1,0],[2,0],[3,0],[4,1],[5,1]], G[[0,0],[1,0],[2,0]] ]
    ;

This output corresponds to the coproduct of the sum of all 1PI,
2-point, 2-loop diagrams in phi^4 theory.

For more exhaustive tests run the scripts in the `tests` directory.

    $ cd tests
    $ .test_feyngen.sh
    $ .test_feyncop.sh

