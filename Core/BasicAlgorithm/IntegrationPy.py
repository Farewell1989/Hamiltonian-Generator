from numpy import *
from numpy.polynomial.legendre import leggauss
def knots_and_weights(a,b,deg,method='legendre'):
    if method=='legendre':
        knots,weights=leggauss(deg)
        knots=(b-a)/2*knots+(a+b)/2
        weights=(b-a)/2*weights
    return knots,weights

def integration(func,a,b,args=(),deg=64,method='legendre'):
    knots,weights=knots_and_weights(a,b,deg,method)
    result=0
    for knot,weight in zip(knots,weights):
        result+=func(knot,*args)*weight
    return result
        
