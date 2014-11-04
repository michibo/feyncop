
""" outputfg.py: This file is part of the feyncop/feyngen package.
    Implements graph output helper functions. """

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


from fractions import *

def get_element_str( element, fac ):
    """Helper function: Adds a factor to a string representation of an object 
        if it is not trivial."""

    sign = "-" if fac < 0 else "+"
    if abs(fac) == Fraction(1,1): 
        return "%s %s" % (sign, element)
    else:
        return "%s %d/%d * %s" % (sign, fac.numerator, fac.denominator, element)


