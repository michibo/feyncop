#!/bin/bash -e

# test_feyngen.sh
#
# This file contains a bash script that tests feyngen
# by comparing its output to precomputed data.

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky, 2024

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

APP_ROOT="$(dirname "$(dirname "$(readlink -fm "$0")")")"

FG=$APP_ROOT/feyngen
WS=$APP_ROOT/weightsum

MAXLOOPS=2
MAXSRCS=4

CMPFILE=$APP_ROOT/tests/verification.out
TESTFILE=$(mktemp)

echo "$TESTFILE"

echo "This file contains a list of numbers produced by feyngen with different parameters" >> $TESTFILE

for degree in "4" "3"
do
    for srcs in `seq 0 $MAXSRCS`
    do
        for loopnum in `seq 0 $MAXLOOPS`
        do
            for cntd in "  " "-c"
            do
                for vtx_2_cntd in "  " "-v"
                do
                    for tadpole in "  " "-t"
                    do
                        for edge_2_cntd in "  " "-p"
                        do
                            for label in "      " "-u"
                            do
                                echo ./feyngen -k$degree $loopnum -j$srcs $cntd $vtx_2_cntd $tadpole $edge_2_cntd $label | tee -a $TESTFILE
                                $FG -k$degree $loopnum -j$srcs $cntd $vtx_2_cntd $tadpole $edge_2_cntd $label | $WS $label | tee -a $TESTFILE
                            done
                        done
                    done
                done
            done
        done
    done
done

MAXBOSONS=2
MAXFERMIONS=2

for bosons in `seq 0 $MAXBOSONS`
do
    for fermions in `seq 0 $MAXFERMIONS`
    do
        for loopnum in `seq 0 $MAXLOOPS`
        do
            for cntd in "  " "-c"
            do
                for vtx_2_cntd in "  " "-v"
                do
                    for tadpole in "  " "-t"
                    do
                        for edge_2_cntd in "  " "-p"
                        do
                            for label in "      " "-u"
                            do
                                echo ./feyngen --qed $loopnum -b$bosons -f$fermions $cntd $vtx_2_cntd $tadpole $edge_2_cntd $label | tee -a $TESTFILE
                                $FG --qed $loopnum -b$bosons -f$fermions $cntd $vtx_2_cntd $tadpole $edge_2_cntd $label | $WS $label | tee -a $TESTFILE
                            done
                        done
                    done
                done
            done
        done
    done
done

MAXLOOPS=1
MAXBOSONS=2
MAXFERMIONS=2
MAXGHOSTS=2

for bosons in `seq 0 $MAXBOSONS`
do
    for fermions in `seq 0 $MAXFERMIONS`
    do
        for ghosts in `seq 0 $MAXGHOSTS`
        do
            for loopnum in `seq 0 $MAXLOOPS`
            do
                for cntd in "  " "-c"
                do
                    for vtx_2_cntd in "  " "-v"
                    do
                        for tadpole in "  " "-t"
                        do
                            for edge_2_cntd in "  " "-p"
                            do
                                for label in "      " "-u"
                                do
                                    echo ./feyngen --ym $loopnum -b$bosons -f$fermions $cntd $vtx_2_cntd $tadpole $edge_2_cntd $label | tee -a $TESTFILE
                                    $FG --ym $loopnum -b$bosons -f$fermions $cntd $vtx_2_cntd $tadpole $edge_2_cntd $label | $WS $label | tee -a $TESTFILE
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done

DIFF=$(diff $TESTFILE $CMPFILE) 

if [ "$DIFF" == "" ] 
then
    echo "TEST ENDED SUCCESSFULLY"
else
    echo "TEST FAILED"
    echo $DIFF
fi

rm $TESTFILE
