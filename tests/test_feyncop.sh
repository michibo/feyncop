#!/bin/bash -e

# test_feyncop.sh
#
# This file contains a bash script that tests feyncop
# by checking a nontrivial identity within the Hopf algebra
# of Feynman graphs

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky, 2024

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email


APP_ROOT="$(dirname "$(dirname "$(readlink -fm "$0")")")"

FG=$APP_ROOT/feyngen
FC=$APP_ROOT/feyncop

$FG -k4 -j4 -pu 2 | $FC -k
$FG -k4 -j4 -pu 3 | $FC -k
$FG -k4 -j4 -pu 4 | $FC -k
$FG -k4 -j4 -pu 5 | $FC -k

$FG -k3 -j3 -pu 2 | $FC -k -D6
$FG -k3 -j3 -pu 3 | $FC -k -D6
$FG -k3 -j3 -pu 4 | $FC -k -D6
