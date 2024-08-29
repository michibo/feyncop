
""" outputfg.py: This file is part of the feyncop/feyngen package.
    Implements graph output helper functions. """

# See also: https://github.com/michibo/feyncop

# Author: Michael Borinsky
# Python3 port: Frédéric Chapoton

# Bugreports, comments, or suggestions are always welcome.
# For instance, via github or email


from fractions import *

def get_element_str( element, fac ):
    """Helper function: Adds a factor to a string representation of an object 
        if it is not trivial."""

    sign = "-" if fac < 0 else "+"
    return "%s %d/%d * %s" % (sign, abs(fac.numerator), fac.denominator, element)


