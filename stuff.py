
"""stuff.py: This file is part of the feyncop/feyngen package.
    Implements some mathematical helper functions. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

from math import comb, factorial


def parse_cntd(cntd, vtx2cntd) -> str:
    """
    Manage what to ask nauty about connectivity.
    """
    if vtx2cntd:
        return "connected"  # should be "biconnected"
    if cntd:
        return "connected"
    return ""


def flip(xy):
    """
    to flip edges

    EXAMPLES::

        sage: from stuff import *
        sage: flip((4, 5))
        (5, 4)
    """
    x, y = xy
    return (y, x)


def double_factorial(k):
    """
    Calculate the double factorial of k.

    EXAMPLES::

        sage: from stuff import *
        sage: double_factorial(-1)
        1
        sage: double_factorial(5)
        15
        sage: double_factorial(6)
        48
    """
    if k == -1:
        return 1

    if k % 2:
        n = (k - 1) // 2
        return factorial(2 * n + 1) // 2**n // factorial(n)

    n = k // 2
    return factorial(n) * 2**n


def binomial(n, k):
    """
    Calculate the binomial coefficient n over k.

    EXAMPLES::

        sage: from stuff import *
        sage: binomial(6,-2)
        0
        sage: binomial(6,2)
        15
        sage: binomial(-6,2)
        21
    """
    if k < 0:
        return 0
    if n < 0:
        return binomial(-n + k - 1, k) * (-1)**k
    return comb(n, k)
