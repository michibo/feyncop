
import collections
from stuff import *
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
#if A[0] != 0:
#           raise

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

def n_from_L( L, k, m ):
    n_prime = (m + 2*(L - 1))
    if n_prime % (k - 2) != 0 : return 0
    return n_prime / (k - 2)

def phi_k_class_coeff( num_loops, num_ext_edges, vtx_degree ):
    num_vtcs = n_from_L( num_loops, vtx_degree, num_ext_edges )

    return phi_k_cc( num_vtcs, num_ext_edges, vtx_degree )

def cntd_phi_k_class_coeff( num_loops, num_ext_edges, vtx_degree ):
    num_vtcs = n_from_L( num_loops, vtx_degree, num_ext_edges )

    return cntd_phi_k_cc( num_vtcs, num_ext_edges, vtx_degree )

def e2c_phi_k_class_coeff( num_loops, num_ext_edges, vtx_degree ):
    num_vtcs = n_from_L( num_loops, vtx_degree, num_ext_edges )

    return e2c_phi_k_cc( num_vtcs, num_ext_edges, vtx_degree )

def qed_class_coeff( num_loops, num_ext_fermions, num_ext_bosons ):
    num_vtcs = n_from_L( num_loops, 3, num_ext_fermions + num_ext_bosons )

    return qed_cc( num_vtcs, num_ext_fermions, num_ext_bosons )

# 1/n! ( g^2 \phi^4 + g \phi^3 + g \phi \psi^2 )^n
# 1/n! \sum_{ k1+k2+k3 = n } ( n \over k1 k2 k3 ) g^( 2*k1 + k2 + k3 ) \phi^( 4*k1 + 3*k2 + k1 ) \psi^(2*k1)
def qcd_class_coeff_no_ghosts( num_loops, num_ext_fermions, num_ext_bosons ):
    if num_ext_fermions % 2 != 0:
        return 0
    L = num_loops
    r = num_ext_fermions/2
    m = num_ext_bosons
    s = m + 2*r + 2*(L - 1)

    return sum( qcd_cc_no_ghosts( l, s, m, r ) for l in range( 0, (3*s + m + 2*r)/2+1 ) )

def qcd_class_coeff( num_loops, num_ext_fermions, num_ext_bosons ):
    if num_ext_fermions % 2 != 0:
        return 0
    L = num_loops
    r = num_ext_fermions/2
    m = num_ext_bosons
    s = m + 2*r + 2*(L - 1)

    return sum( qcd_cc( l, s, m, r ) for l in range( 0, (3*s + m + 2*r)/2+1 ) )


def qcd_cc_no_ghosts( l, s, m, r ):
    max_q = max_p = (2*l - m - 2*r - 2*s) 

    def gen_terms():
        for q in range( 0, max_q+1 ):
            for p in range( 0, 1 ):
                n = 2*l - 2*s - q - p - m - 2*r
                two_k = 3*s + m + 2*r - 2*l
                if (two_k < 0) or (n < 0) or (two_k % 2 != 0):
                    continue
                k = two_k / 2
                
                phi_power = 3*n + 4*k + q + p + m 
                psi_power = 2*r + 2*q
                c_power = 2*p

                exp_coeff = factorial(n) * factorial(k) * factorial(q) * \
                            factorial(p) * factorial(m) * factorial(r)**2
                dom1_coeff = factorial(3) ** n
                dom2_coeff = factorial(4) ** k
                wick = double_factorial( phi_power - 1 ) * factorial( psi_power / 2) * factorial( c_power / 2)
                
                yield Fraction( wick, exp_coeff*dom1_coeff*dom2_coeff )

    return sum( fac for fac in gen_terms() )

def qcd_cc( l, s, m, r ):
    max_q = max_p = (2*l - m - 2*r - 2*s) 

    def gen_terms():
        for q in range( 0, max_q+1 ):
            for p in range( 0, max_p+1 ):
                n = 2*l - 2*s - q - p - m - 2*r
                two_k = 3*s + m + 2*r - 2*l
                if (two_k < 0) or (n < 0) or (two_k % 2 != 0):
                    continue
                k = two_k / 2
                
                phi_power = 3*n + 4*k + q + p + m 
                psi_power = 2*r + 2*q
                c_power = 2*p

                exp_coeff = factorial(n) * factorial(k) * factorial(q) * \
                            factorial(p) * factorial(m) * factorial(r)**2
                dom1_coeff = factorial(3) ** n
                dom2_coeff = factorial(4) ** k
                wick = double_factorial( phi_power - 1 ) * factorial( psi_power / 2) * factorial( c_power / 2)
                
                yield Fraction( wick, exp_coeff*dom1_coeff*dom2_coeff )

    return sum( fac for fac in gen_terms() )

def cntd_qed_class_coeff( num_loops, num_ext_fermions, num_ext_bosons):
    num_vtcs = n_from_L( num_loops, 3, num_ext_fermions + num_ext_bosons )

    return cntd_qed_cc( num_vtcs, num_ext_fermions, num_ext_bosons )

def cntd_qcd_class_coeff_no_ghosts( num_loops, num_ext_fermions, num_ext_bosons ):
    if num_ext_fermions % 2 != 0:
        return 0
    l = num_loops
    r = num_ext_fermions/2
    m = num_ext_bosons
    s = m + 2*r + 2*(l - 1)

    return cntd_qcd_cc_no_ghosts( s, m, r )

def cntd_qcd_class_coeff( num_loops, num_ext_fermions, num_ext_bosons ):
    if num_ext_fermions % 2 != 0:
        return 0
    l = num_loops
    r = num_ext_fermions/2
    m = num_ext_bosons
    s = m + 2*r + 2*(l - 1)

    return cntd_qcd_cc( s, m, r )

def phi_k_cc( n, m, k ):
    dbl_l = m + n * k

    if dbl_l % 2 != 0:
        return 0
    
    exp_coeff = factorial(n) * factorial(m)
    wick_coeff = double_factorial( dbl_l - 1 )

    lambda_coeff = factorial( k ) ** n

    return Fraction( wick_coeff, lambda_coeff*exp_coeff )

def phi_3_4_cc( n, m, k ):
    dbl_l = m + n * k

    if dbl_l % 2 != 0:
        return 0
    
    exp_coeff = factorial(n) * factorial(m)
    wick_coeff = double_factorial( dbl_l - 1 )

    lambda_coeff = factorial( k ) ** n

    return Fraction( wick_coeff, lambda_coeff*exp_coeff )


def phi_k_cc_l( l, m, k ):
    return phi_k_cc( (2*l - m)/k, m, k ) if (2*l - m) % k == 0 else 0 

def cntd_phi_k_cc( n, m, k ):
    dbl_l = m + n * k

    if dbl_l % 2 != 0:
        return 0

    l = dbl_l/2

    A = [ [ phi_k_cc_l(o, p, k) for p in range(2*o + 1) ] for o in range(l+1) ]
    Alog = lLog( A )
    
    return Alog[l][m]

def e2c_phi_k_cc( n, m, k ):
    dbl_l = m + n * k

    if dbl_l % 2 != 0:
        return 0

    l = dbl_l/2

    A = [ [ phi_k_cc_l(o, p, k) if p <= 2*o else 0 for p in range(2*l + 1) ] for o in range(l+1) ]
    Alog = lLog( A )
    
    W = [ [ Alog[o][p] for o in range(l+1) ] for p in range(2*l + 1) ]

    phi_c = lDerivative( W )

    j = lCompInverse( phi_c )
    
    Gamma = [0] + lSum( lComposition( W, j )[1:], lScalMult(-1, phi_c) )

    return Gamma[m][l]

def qed_cc( num_vtcs, num_ext_fermions, num_ext_bosons):
    if num_ext_fermions % 2 != 0:
        return 0
    phi_power = num_vtcs + num_ext_bosons
    psi_power = num_ext_fermions + 2*num_vtcs

    if phi_power % 2 != 0:
        return 0

    exp_coeff = factorial( num_vtcs ) * factorial(num_ext_bosons) * factorial( num_ext_fermions/2 )**2
    wick_coeff1 = double_factorial( phi_power - 1 )
    wick_coeff2 = factorial( psi_power/2 )

    return Fraction( wick_coeff1*wick_coeff2, exp_coeff )

def qed_cc_l( l, m, r ):
    return qed_cc( (2*l - m - r)/3, m, r ) if (2*l - m - r) % 3 == 0 else 0

def cntd_qed_cc( num_vtcs, num_ext_bosons, num_ext_fermions ):
    num_edges_dbl = 3*num_vtcs + num_ext_bosons + num_ext_fermions
    if num_edges_dbl % 2 != 0:
        return 0

    num_edges = num_edges_dbl / 2
    A = [ [ [ qed_cc_l(l, m, r) for m in range(2*l - r + 1) ] for r in range(2*l+1) ] for l in range(num_edges+1) ]
    Alog = lLog( A )
    
    return Alog[num_edges][num_ext_fermions][num_ext_bosons]

def cntd_qcd_cc_no_ghosts( s, m, r ):
    l = (3*s + m + 2*r)/2
    # Warning: quick and dirty optimization!
    A = [ [ [ [ qcd_cc_no_ghosts( lp, sp, mp, rp ) for sp in range(max( min(s+1, lp+1-(mp+2*rp)/2), 1)) ] for mp in range(m+1) ] for rp in range(r+1) ] for lp in range(l+1) ]

    Alog = lLog( A )
    return sum( Alog[lp][r][m][s] for lp in range(l+1) if len(Alog[lp][r][m]) > s )

def cntd_qcd_cc( s, m, r ):
    l = (3*s + m + 2*r)/2
    # Warning: quick and dirty optimization!
    A = [ [ [ [ qcd_cc( lp, sp, mp, rp ) for sp in range(max( min(s+1, lp+1-(mp+2*rp)/2), 1)) ] for mp in range(m+1) ] for rp in range(r+1) ] for lp in range(l+1) ]

    Alog = lLog( A )
    return sum( Alog[lp][r][m][s] for lp in range(l+1) if len(Alog[lp][r][m]) > s )
