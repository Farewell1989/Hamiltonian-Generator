from numpy import *
from scipy.linalg import eigh
def berry_curvature(H,kx,ky,mu,d=10**-6):
    '''
    This function calculates the Berry curature of the occupied bands for a Hamiltonian with the given chemical potential.
    '''
    result=0
    Vx=(H(kx+d,ky)-H(kx-d,ky))/(2*d)
    Vy=(H(kx,ky+d)-H(kx,ky-d))/(2*d)
    Es,Evs=eigh(H(kx,ky))
    for n in xrange(Es.shape[0]):
        for m in xrange(Es.shape[0]):
            if Es[n]<=mu and Es[m]>mu:
                result-=2*(vdot(dot(Vx,Evs[:,n]),Evs[:,m])*vdot(Evs[:,m],dot(Vy,Evs[:,n]))/(Es[n]-Es[m])**2).imag
    return result
