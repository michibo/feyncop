
"""powerseries.py: This file is part of the feyncop/feyngen package.
    Collection of subroutines for the manipulation of multivariable polynomials, which can be seen as truncated multivariable power series."""

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

from fractions import Fraction
from functools import reduce


def unary_rec_list_op( op, A ):
    """Apply an unary operation to a multivariable polynomial: op(A)"""

    if type(A) is list:
        return [ unary_rec_list_op( op, a ) for a in A ]
    else:
        return op( A )

def binary_rec_list_op( op, A, B ):
    """Apply a binary operation to two multivariable polynomials: op(A,B)"""

    if type(A) is list and type(B) is list:
        return [ binary_rec_list_op( op, a, b ) for a,b in zip(A,B) ]
    else:
        return op( A, B )

def lSum( A, B ):
    """Sum two multivariable polynomials: A+B"""

    def help_sum( a, b ): return a+b

    return binary_rec_list_op( help_sum, A, B )

def lScalMult( m, A ):
    """Scalar multiply a multivariable polynomials: m*A"""

    def help_mul( a ): return m*a

    return unary_rec_list_op( help_mul, A )

def lConvolute( A, B ):
    """Multiply/Convolute two multivariable polynomials: A*B"""

    if type(A) is list and type(B) is list:
        return [ reduce( lSum, ( lConvolute(A[k], B[n-k]) for k in range(n+1) if k < len(A) and (n-k) < len(B)) ) for n in range( len(A) + len(B) - 1) ]
    else:
        return A*B

def lInvert( A ):
    """Calculate reciproke truncated power series: 1/A"""

    if type(A) is list:
        if len(A) > 1:
            Ainv_s = lInvert( A[:-1] )
            Ap = [ reduce( lSum, ( lConvolute( Ainv_s[k], A[n-k] ) for k in range(n) ) ) for n in range(1, len(A) ) ]
            A0rec = lInvert(A[0])
            A0rec_neg = lScalMult( -1, A0rec )
            return [ A0rec ] + [ lConvolute( A0rec_neg, a ) for a in Ap ]
        else:
            return [ lInvert(A[0]) ]
    else:
        return Fraction(1, A)

def lLog( A ):
    """Calculate the log of A: log(A)"""

    if type(A) is list:
        Ainv = lInvert( A )
        Ap = [ lScalMult(Fraction(1, n), reduce( lSum, ( lScalMult( k, lConvolute( A[k], Ainv[n-k] ) ) for k in range(1,n+1) ) ) ) for n in range(1,len(A)) ]
        return [lLog(A[0])] + Ap
    else:
        if A == 1:
            return 0
        else:
            return log( A )
