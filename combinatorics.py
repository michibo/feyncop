
"""combinatorics.py: This file is part of the feyncop/feyngen package.
    Collection of subroutines for the combinatorial calculation of zero dimensional field theories."""

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

from stuff import *
from fractions import Fraction
from powerseries import lLog,lConvolute

def phi_k_class_coeff( L, m, k ):
    """Calculate the sum of the symmetry factors of all phi^k diagrams 
        with L loops, m external edges and valency k."""
    
    s_tkm2 = m + 2*(L - 1)
    if s_tkm2 % (k-2) != 0: return 0 
    s = s_tkm2//(k-2)

    if s<0:
        return 0

    return phi_k_cc( s, m, k )

def phi34_class_coeff( L, m ):
    """Calculate the sum of the symmetry factors of all (phi^3 + phi^4)
        diagrams with L loops, m external edges."""

    s = 2*(L-1) + m

    if s<0:
        return 0

    return phi34_cc( s, m )

def qed_class_coeff( L, rt2, m ):
    """Calculate the sum of the symmetry factors of all QED diagrams 
        with L loops, rt2 external fermions and m external bosons."""

    if rt2 % 2 != 0:
        return 0
    
    r = rt2//2
    s = m + 2*r + 2*(L - 1)
    
    if s<0:
        return 0

    return qed_cc(s,m,r)

def qed_furry_class_coeff( L, rt2, m ):
    """Calculate the sum of the symmetry factors of all QED diagrams 
        respecting Furry's theorem
        with L loops, rt2 external fermions and m external bosons."""

    if rt2 % 2 != 0:
        return 0
    
    r = rt2//2
    s = m + 2*r + 2*(L - 1)
    
    if s<0:
        return 0

    return qed_furry_cc(s,m,r)

def qcd_class_coeff( L, rt2, ut2, m ): 
    """Calculate the sum of the symmetry factors of all QCD diagrams 
       with L loops, rt2 external fermions, ut2 external ghosts and 
       m external bosons."""

    if rt2 % 2 != 0:
        return 0
    if ut2 % 2 != 0:
        return 0
    
    r = rt2//2
    u = ut2//2
    s = m + 2*r + 2*u + 2*(L - 1)

    if s<0:
        return 0

    return qcd_cc(s,m,r,u)

def cntd_phi_k_class_coeff( L, m, k ):
    """Calculate the sum of the symmetry factors of all connected phi^k 
        diagrams with L loops, m external edges and valency k."""
    
    s_tkm2 = m + 2*(L - 1)
    if s_tkm2 % (k-2) != 0: return 0 
    s = s_tkm2//(k-2)

    if s<0:
        return 0

    A = [ [ phi_k_cc(sp, mp, k) for mp in range(m+1) ] for sp in range(s+1) ]
    Alog = lLog( A )
    
    return Alog[s][m]

def cntd_phi34_class_coeff( L, m ):
    """Calculate the sum of the symmetry factors of all connected 
        (phi^3 + phi^4) diagrams with L loops, m external edges."""

    s = m + 2*(L - 1)

    if s<0:
        return 0

    A = [ [ phi34_cc( sp, mp ) for mp in range(m+1) ] for sp in range(s+1) ]

    Alog = lLog( A )
    return Alog[s][m]

def cntd_qed_class_coeff( L, rt2, m):
    """Calculate the sum of the symmetry factors of all connected 
        QED diagrams with L loops, rt2 external fermions and 
        m external bosons."""

    if rt2 % 2 != 0:
        return 0
    
    r = rt2//2
    s = m + 2*r + 2*(L - 1)

    if s<0:
        return 0

    A = [ [ [ qed_cc(sp, mp, rp) for mp in range(m+1) ] for rp in range(r+1) ] for sp in range(s+1) ]
    Alog = lLog( A )
    
    return Alog[s][r][m]

def cntd_qed_furry_class_coeff( L, rt2, m):
    """Calculate the sum of the symmetry factors of all connected 
        QED diagrams respecting Furry's theorem
        with L loops, rt2 external fermions and 
        m external bosons."""

    if rt2 % 2 != 0:
        return 0
    
    r = rt2//2
    s = m + 2*r + 2*(L - 1)

    if s<0:
        return 0

    A = [ [ [ qed_furry_cc(sp, mp, rp) for mp in range(m+1) ] for rp in range(r+1) ] for sp in range(s+1) ]
    Alog = lLog( A )
    
    return Alog[s][r][m]
def cntd_qcd_class_coeff( L, rt2, ut2, m ):
    """Calculate the sum of the symmetry factors of all connected QCD 
        diagrams with L loops, rt2 external fermions, ut2 external ghosts 
        and m external bosons."""

    if rt2 % 2 != 0:
        return 0
    
    r = rt2//2
    u = ut2//2
    s = m + 2*r + 2*u + 2*(L - 1)

    if s<0:
        return 0

    A = [ [ [ [ qcd_cc( sp, mp, rp, up ) for rp in range(r+1) ] for up in range(u+1) ] for mp in range(m+1) ] for sp in range(s+1) ] 

    Alog = lLog( A )
    return Alog[s][m][u][r]

def phi_k_cc( s, m, k ):
    """Helper function which evaluates the relevant term in the 
        generating function(al) of a zero dimensional phi^k theory. 
        s corresponds to the power in the coupling constant and 
        m to the power of the field sources."""
        
    l_t2 = m + s * k

    if l_t2 % 2 != 0:
        return 0
    
    denom = factorial(s) * factorial(m) * factorial( k ) ** s
    nom = double_factorial( l_t2 - 1 )

    return Fraction( nom, denom )

def phi34_cc( s, m ):
    """Helper function which evaluates the relevant term in the 
        generating function(al) of a zero dimensional (phi^3+phi^4) theory. 
        s corresponds to the power in the coupling constant and 
        m to the power of the field sources."""

    l_min = max(0, (m + 2*s + 1) // 2)
    l_max = ( 3*s + m ) // 2 
    S = 0
    for l in range(l_min, l_max+1):
        n_1 =  2*l - 2*s - m
        n_2_t2 = -2*l + 3*s + m
        if n_2_t2 % 2 != 0:
            continue
        n_2 = n_2_t2 // 2

        denom = factorial( n_1 ) * factorial( n_2 ) * \
                factorial( m ) * factorial(3)**n_1 *  \
                factorial(4) ** n_2

        S += Fraction(double_factorial(2*l-1), denom)

    return S

def qed_cc( s, m, r):
    """Helper function which evaluates the relevant term in the 
        generating function(al) of zero dimensional QED. 
        s corresponds to the power in the coupling constant, 
        m to the power of the photon field sources and r to 
        the power of the absolute squared fermion sources."""

    l1_t2 = s + m
    l2 = r + s

    if l1_t2 % 2 != 0:
        return 0

    denom = factorial(s) * factorial(m) * factorial(r)**2
    nom = double_factorial( l1_t2 - 1 ) * factorial( l2 )

    return Fraction( nom, denom )

def qed_furry_cc( s, m, r):
    """Helper function which evaluates the relevant term in the 
        generating function(al) of zero dimensional QED respecting 
        Furry's theorem. 
        s corresponds to the power in the coupling constant, 
        m to the power of the photon field sources and r to 
        the power of the absolute squared fermion sources."""
    
    S1 = [ binomial( r + n - 1, n ) for n in range(s+1) ]
    S2 = [ Fraction(binomial(2*n,n), 4**n) for n in range(s+1) ]
    S3 = [ (-1)**n * Fraction(binomial(2*n,n), 4**n) for n in range(s+1) ]

    C = lConvolute(lConvolute( S1, S2 ), S3)

    l1_t2 = s + m
    l2 = r + s

    if l1_t2 % 2 != 0:
        return 0

    denom = factorial(m) * factorial(r)
    nom = double_factorial( l1_t2 - 1 )

    return Fraction( nom * C[s], denom )

def qcd_cc( s, m, r, u ):
    """Helper function which evaluates the relevant term in the 
        generating function(al) of zero dimensional QCD. 
        s corresponds to the power in the coupling constant, 
        m to the power of the photon field sources and r/u to 
        the power of the absolute squared fermion/ghost sources."""

    l2_min = r
    l3_min = u
    l1_min = (m+1)//2
    l1_max = l2_max = l3_max = (3*s+m+2*r+2*u)//2

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
                n2 = n2_t2//2
                if n1 < 0 or n2 < 0 or n3 < 0 or n4 < 0:
                    continue

                denom = factorial(n1)*factorial(n2)*factorial(n3)*factorial(n4)*factorial(3)**n1*factorial(4)**n2*factorial(m)*factorial(r)**2*factorial(u)**2

                nom = double_factorial(2*l1-1)*factorial(l2)*factorial(l3)
                S+= Fraction(nom, denom)

    return S
