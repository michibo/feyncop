
"""stuff.py: Implements some mathematical helper functions. """

__author__ = "Michael Borinsky"
__email__ = "borinsky@physik.hu-berlin.de"



from math import *

def double_factorial( k ):
    """Calculate the double factorial of k."""
    if k == -1:
        return 1
    if k % 2 == 1:
        n = (k - 1) / 2 
        return factorial( 2*n + 1 ) / 2**n / factorial(n)
    else:
        n = k / 2
        return factorial( n ) * 2**n

def binomial( n, k ):
    """Calculate the binomial coefficient n over k."""
    if k < 0:
        return 0
    if n < 0:
        return binomial( -n + k - 1, k ) * (-1)**k

    b = 1
    for i in range(k):
        b *= (n-i)
        b /= (1+i)
    return b

