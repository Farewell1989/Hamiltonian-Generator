from numpy import *
from numpy.polynomial.legendre import leggauss
def integration_knots_weights(a,b,deg,method='legendre'):
    if method=='legendre':
        knots,weights=leggauss(deg)
        knots=(b-a)/2*knots+(a+b)/2
        weights=(b-a)/2*weights
    return knots,weights
