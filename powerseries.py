
"""powerseries.py: This file is part of the feyncop/feyngen package.
    Collection of subroutines for the manipulation of multivariable polynomials, which can be seen as truncated multivariable power series."""

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email

from fractions import Fraction
from functools import reduce
from math import log


def unary_rec_list_op(op, A):
    """Apply an unary operation to a multivariable polynomial: op(A)"""

    if type(A) is list:
        return [unary_rec_list_op(op, a) for a in A]
    else:
        return op(A)


def binary_rec_list_op(op, A, B):
    """Apply a binary operation to two multivariable polynomials: op(A,B)"""

    if isinstance(A, list) and isinstance(B, list):
        return [binary_rec_list_op(op, a, b) for a, b in zip(A, B)]
    return op(A, B)


def lSum(A, B):
    """Sum two multivariable polynomials: A+B"""

    def help_sum(a, b):
        return a + b

    return binary_rec_list_op(help_sum, A, B)


def lScalMult(m, A):
    """Scalar multiply a multivariable polynomials: m*A"""

    def help_mul(a):
        return m * a

    return unary_rec_list_op(help_mul, A)


def lConvolute(A, B):
    """Multiply/Convolute two multivariable polynomials: A*B"""

    if isinstance(A, list) and isinstance(B, list):
        return [reduce(lSum, (lConvolute(A[k], B[n - k]) for k in range(n + 1)
                              if k < len(A) and (n - k) < len(B)))
                for n in range(len(A) + len(B) - 1)]
    return A * B


def lInvert(A):
    """Calculate reciproke truncated power series: 1/A"""

    if type(A) is list:
        if len(A) > 1:
            Ainv_s = lInvert(A[:-1])
            Ap = [reduce(lSum, (lConvolute(Ainv_s[k], A[n - k]) for k in range(n))) for n in range(1, len(A))]
            A0rec = lInvert(A[0])
            A0rec_neg = lScalMult(-1, A0rec)
            return [A0rec] + [lConvolute(A0rec_neg, a) for a in Ap]
        else:
            return [lInvert(A[0])]
    else:
        return Fraction(1, A)


def lLog(A):
    """Calculate the log of A: log(A)"""

    if isinstance(A, list):
        Ainv = lInvert(A)
        Ap = [lScalMult(Fraction(1, n),
                        reduce(lSum, (lScalMult(k, lConvolute(A[k], Ainv[n - k]))
                                      for k in range(1, n + 1))))
              for n in range(1, len(A))]
        return [lLog(A[0])] + Ap

    if A == 1:
        return 0

    return log(A)
