
from math import *

def wick_contractions( l ):
    if len(l) == 1:
        return
    if len(l) == 0:
        yield []
        return

    for i in range(1, len(l)):
        for ctrctn in wick_contractions( l[1:i] + l[i+1:] ):
            yield [(l[0],l[i])] + ctrctn
        

def double_factorial( k ):
    if k == -1:
        return 1
    if k % 2 == 1:
        n = (k - 1) / 2 
        return factorial( 2*n + 1 ) / 2**n / factorial(n)
    else:
        n = k / 2
        return factorial( n ) * 2**n

def binomial( n, k ):
    if k < 0:
        return 0
    if n < 0:
        return binomial( -n + k - 1, k ) * (-1)**k

    b = 1
    for i in range(k):
        b *= (n-i)
        b /= (1+i)
    return b

def multinomial( ks ):
    nom = sum(ks)

    res = factorial( nom )

    for k in ks:
        res/= factorial(k)

    return res

def ordered_partitions_with_0( n, m ):
    for partition in partitions(n):
        if len(partition) > m:
            continue
        partition += [0]*(len(partition)-m)

        for perm in perm_unique( partition ):
            yield perm

def partitions(n):
    if n == 0:
        yield ()
        return

    for p in partitions(n-1):
        yield (1,) + p
        if p and (len(p) < 2 or p[1] > p[0]):
            yield (p[0] + 1,) + p[1:]

class unique_element:
    def __init__(self,value,occurrences):
        self.value = value
        self.occurrences = occurrences

def perm_unique(elements):
    eset=set(elements)
    listunique = [unique_element(i,elements.count(i)) for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)

def perm_unique_helper(listunique,result_list,d):
    if d < 0:
        yield tuple(result_list)
    else:
        for i in listunique:
            if i.occurrences > 0:
                result_list[d]=i.value
                i.occurrences-=1
                for g in  perm_unique_helper(listunique,result_list,d-1):
                    yield g
                i.occurrences+=1

