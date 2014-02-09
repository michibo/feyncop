
from fractions import Fraction

def unary_rec_list_op( op, A ):
    if type(A) is list:
        return [ unary_rec_list_op( op, a ) for a in A ]
    else:
        return op( A )

def binary_rec_list_op( op, A, B ):
    if type(A) is list and type(B) is list:
        return [ binary_rec_list_op( op, a, b ) for a,b in zip(A,B) ]
    else:
        return op( A, B )

def lSum( A, B ):
    def help_sum( a, b ): return a+b

    return binary_rec_list_op( help_sum, A, B )

def lScalMult( m, A ):
    def help_mul( a ): return m*a

    return unary_rec_list_op( help_mul, A )

def lConvolute( A, B ):
    if type(A) is list and type(B) is list:
        return [ reduce( lSum, ( lConvolute(A[k], B[n-k]) for k in range(n+1) if k < len(A) and (n-k) < len(B)) ) for n in range( len(A) + len(B) - 1) ]
    else:
        return A*B

def lComposition( A, B ):
    if type(A) is list:

        def gen_coeffs():
            for n in range(1, len(A)):
                def tens_prod_op( b ):
                    return lScalMult( b, A[n] )
                yield unary_rec_list_op( tens_prod_op, lPower(B, n) )

        return reduce( lSum, gen_coeffs() )
    else:
        raise

def lCompInverse( A ):
    if A[0] != 0:
        raise

    Phi = lInvert(A[1:])

    return [0] + [ lScalMult(Fraction(1,n), lPower(Phi, n)[n-1]) for n in range(1, len(A)) ]

def lDerivative( A ):
    if type(A) is list:
        return [ lScalMult( n, A[n-1] ) for n in range(1,len(A)+1) ]
    else:
        return 0

def lPower( A, n ):
    if type(A) is list:
        def gen_pow( m ):
            for partition in partitions( m ):
                cntr = collections.Counter( partition )
                p = [ cntr[i] for i in range(1,m+1) ]
                
                if sum(p) > n : continue

                p = [ n - sum(p) ] + p

                fac = multinomial( p )

                yield lScalMult( fac, reduce( lConvolute, ( lPower( a, l ) for a,l in zip(A, p) if l > 0 ) ) )

        return [ reduce( lSum, gen_pow(m) ) for m in range( len(A) ) ]
    else:
        return A ** n

def lInvert( A ):
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
    if type(A) is list:
        Ainv = lInvert( A )
        Ap = [ lScalMult(Fraction(1, n), reduce( lSum, ( lScalMult( k, lConvolute( A[k], Ainv[n-k] ) ) for k in range(1,n+1) ) ) ) for n in range(1,len(A)) ]
        return [lLog(A[0])] + Ap
    else:
        if A == 1:
            return 0
        else:
            return log( A )
