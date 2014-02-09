
import collections
from stuff import *
from fractions import Fraction
from powerseries import lLog

def phi_k_class_coeff( num_loops, num_ext_edges, vtx_degree ):
    L = num_loops
    m = num_ext_edges
    k = vtx_degree

    s_tkm2 = m + 2*(L - 1)
    if s_tkm2 % (k-2) != 0: return 0 
    s = s_tkm2/(k-2)

    return phi_k_cc( s, m, k )

def cntd_phi_k_class_coeff( num_loops, num_ext_edges, vtx_degree ):
    L = num_loops
    m = num_ext_edges
    k = vtx_degree

    s_tkm2 = m + 2*(L - 1)
    if s_tkm2 % (k-2) != 0: return 0 
    s = s_tkm2/(k-2)

    return cntd_phi_k_cc( s, m, k )
   
def phi34_class_coeff( num_loops, num_ext_edges ):
    L = num_loops
    m = num_ext_edges

    s = 2*(L-1) + m

    return phi34_cc( s, m )

def qed_class_coeff( num_loops, num_ext_fermions, num_ext_bosons ):
    if num_ext_fermions % 2 != 0:
        return 0
    L = num_loops
    r = num_ext_fermions/2
    m = num_ext_bosons
    s = m + 2*r + 2*(L - 1)

    return qed_cc(s,m,r)

def qcd_class_coeff( num_loops, num_ext_fermions, num_ext_ghosts, num_ext_bosons ):
    if num_ext_fermions % 2 != 0:
        return 0
    if num_ext_ghosts % 2 != 0:
        return 0
    L = num_loops
    r = num_ext_fermions/2
    u = num_ext_ghosts/2
    m = num_ext_bosons
    s = m + 2*r + 2*u + 2*(L - 1)

    return qcd_cc(s,m,r,u)

def phi_k_cc( s, m, k ):
    l_t2 = m + s * k

    if l_t2 % 2 != 0:
        return 0
    
    denom = factorial(s) * factorial(m) * factorial( k ) ** s
    nom = double_factorial( l_t2- 1 )

    return Fraction( nom, denom )

def phi34_cc( s, m ):
    l_min = max(0, (m + 2*s + 1)/2)
    l_max = ( 3*s + m ) / 2 
    S = 0
    for l in range(l_min, l_max+1):
        n_1 =  2*l - 2*s - m
        n_2_t2 = -2*l + 3*s + m
        if n_2_t2 % 2 != 0:
            continue
        n_2 = n_2_t2 / 2

        denom = factorial( n_1 ) * factorial( n_2 ) * \
                factorial( m ) * factorial(3)**n_1 *  \
                factorial(4) ** n_2

        S += Fraction(double_factorial(2*l-1), denom)

    return S

def qed_cc( s, m, r):
    l1_t2 = s + m
    l2 = r + s

    if l1_t2 % 2 != 0:
        return 0

    denom = factorial(s) * factorial(m) * factorial(r)**2
    nom = double_factorial( l1_t2 - 1 ) * factorial( l2 )

    return Fraction( nom, denom )

def qcd_cc( s, m, r, u ):
    l2_min = r
    l3_min = u
    l1_min = (m+1)/2
    l1_max = l2_max = l3_max = (3*s+m+2*r+2*u)/2

    S = 0
    for l1 in range(l1_min, l2_max+1):
        for l2 in range(l2_min, l2_max+1):
            for l3 in range(l3_min, l3_max+1):
                n1 = 2*l1 + l2 + l3 - 2*s - m - r - u
                n2_t2 = -2*(l1+l2+l3) + 3*s + m + 2*r + 2*u
                n3 = l2-r
                n4 = l3-u
                if n2_t2%2 != 0:
                    continue
                n2 = n2_t2/2
                if n1 < 0 or n2 < 0 or n3 < 0 or n4 < 0:
                    continue

                denom = factorial(n1)*factorial(n2)*factorial(n3)*factorial(n4)*factorial(3)**n1*factorial(4)**n2*factorial(m)*factorial(r)**2*factorial(u)**2

                nom = double_factorial(2*l1-1)*factorial(l2)*factorial(l3)
                S+= Fraction(nom, denom)

    return S
                
def cntd_qed_class_coeff( num_loops, num_ext_fermions, num_ext_bosons):
    if num_ext_fermions % 2 != 0:
        return 0
    L = num_loops
    r = num_ext_fermions/2
    m = num_ext_bosons
    s = m + 2*r + 2*(L - 1)

    return cntd_qed_cc(s,m,r)

def cntd_qcd_class_coeff( num_loops, num_ext_fermions, num_ext_ghosts, num_ext_bosons ):
    if num_ext_fermions % 2 != 0:
        return 0
    l = num_loops
    r = num_ext_fermions/2
    u = num_ext_ghosts/2
    m = num_ext_bosons
    s = m + 2*r + 2*u + 2*(l - 1)

    return cntd_qcd_cc( s, m, r, u )

def cntd_phi34_class_coeff( num_loops, num_ext_bosons ):
    L = num_loops
    m = num_ext_bosons
    s = m + 2*(L - 1)

    return cntd_phi34_cc( s, m )

def cntd_phi_k_cc( s, m, k ):
    A = [ [ phi_k_cc(sp, mp, k) for mp in range(m+1) ] for sp in range(s+1) ]
    Alog = lLog( A )
    
    return Alog[s][m]

def cntd_phi34_cc( s, m ):
    A = [ [ phi34_cc( sp, mp ) for mp in range(m+1) ] for sp in range(s+1) ]

    Alog = lLog( A )
    return Alog[s][m]

def cntd_qed_cc( s, m, r ):
    A = [ [ [ qed_cc(sp, mp, rp) for mp in range(m+1) ] for rp in range(r+1) ] for sp in range(s+1) ]
    Alog = lLog( A )
    
    return Alog[s][r][m]

def cntd_qcd_cc( s, m, r, u ):
    A = [ [ [ [ qcd_cc( sp, mp, rp, up ) for rp in range(r+1) ] for up in range(u+1) ] for mp in range(m+1) ] for sp in range(s+1) ] 

    Alog = lLog( A )
    return Alog[s][m][u][r]
