from numpy import *
from scipy.linalg import eigh
def berry_curvature(H,Vx,Vy,mu):
    '''
    This function calculates the Berry curvature of the occupied bands.
    '''
    result=0
    Es,Evs=eigh(H)
    for n in xrange(Es.shape[0]):
        for m in xrange(Es.shape[0]):
            if Es[n]<=mu and Es[m]>mu:
                result-=2*(vdot(dot(Vx,Evs[:,n]),Evs[:,m])*vdot(Evs[:,m],dot(Vy,Evs[:,n]))/(Es[n]-Es[m])**2).imag
    return result
