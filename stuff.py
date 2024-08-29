
"""stuff.py: This file is part of the feyncop/feyngen package.
    Implements some mathematical helper functions. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

from math import *


def double_factorial(k):
    """Calculate the double factorial of k."""
    if k == -1:
        return 1
    if k % 2 == 1:
        n = (k - 1) // 2
        return factorial(2*n + 1) // 2**n // factorial(n)
    else:
        n = k // 2
        return factorial(n) * 2**n


def binomial(n, k):
    """Calculate the binomial coefficient n over k."""
    if k < 0:
        return 0
    if n < 0:
        return binomial(-n + k - 1, k) * (-1)**k

    b = 1
    for i in range(k):
        b *= (n-i)
        b //= (1+i)
    return b
