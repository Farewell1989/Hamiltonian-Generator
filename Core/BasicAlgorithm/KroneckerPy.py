# -*- coding: utf-8 -*-
"""
Created on Sat Jan 30 20:41:33 2016

@author: DZY
"""

from numpy import kron
def kronecker(*arg):
    n=len(arg)
    out=arg[-1]
    for x in xrange(-2,-n-1,-1):
        out=kron(arg[x],out)
    return out