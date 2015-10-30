from numpy import *
from numpy.linalg import norm
from GlobalPy import RZERO

def azimuthd(self):
    '''
    Azimuth in degrees of an array-like vector.
    '''
    if self[1]>=0:
        return degrees(arccos(self[0]/norm(self)))
    else:
        return 360-degrees(arccos(self[0]/norm(self)))

def azimuth(self):
    '''
    Azimuth in radians of an array-like vector.
    '''
    if self[1]>=0:
        return arccos(self[0]/norm(self))
    else:
        return 2*pi-arccos(self[0]/norm(self))

def polard(self):
    '''
    Polar angle in degrees of an array-like vector.
    '''
    if self.shape[0]==3:
        return degrees(arccos(self[2]/norm(self)))
    else:
        raise ValueError("PolarD error: the array-like vector must contain three elements.")

def polar(self):
    '''
    Polar angle in radians of an array-like vector.
    '''
    if self.shape[0]==3:
        return arccos(self[2]/norm(self))
    else:
        raise ValueError("Polar error: the array-like vector must contain three elements.")

def volume(O1,O2,O3):
    '''
    Volume spanned by three array-like vectors.
    '''
    if O1.shape[0] in [1,2] or O2.shape[0] in [1,2] or O3.shape[0] in [1,2]:
        return 0
    elif O1.shape[0] ==3 and O2.shape[0]==3 and O3.shape[0]==3:
        return inner(O1,cross(O2,O3))
    else:
        raise ValueError("Volume error: the shape of the array-like vectors is not supported.")

def is_parallel(O1,O2):
    '''
    Judge whether two array-like vectors are parallel to each other: '0' means not parallel, '1' means parallel, and '-1' means anti-parallel.
    '''
    norm1=norm(O1)
    norm2=norm(O2)
    if norm1<RZERO or norm2<RZERO:
        return 1
    elif O1.shape[0]==O2.shape[0]:
        buff=inner(O1,O2)/(norm1*norm2)
        if abs(buff-1)<RZERO:
            return 1
        elif abs(buff+1)<RZERO:
            return -1
        else:
            return 0
    else:
        raise ValueError("Is_parallel error: the shape of the array-like vectors does not match.") 

class Point:
    '''
    Structured point, which has the following attribute:
    1) site: site index of the point, start with 0;
    2) rcoord: the coordinate in real space;
    3) icoord: the coordinate in lattice space, default value array([]);
    4) atom: atom species on this point, default value 0;
    5) norbital: number of orbitals, default value 1;
    6) nspin: number of spins, default value 2;
    7) nnambu: '1' means not using Nambu space while '2' means using Nambu space, default value 1;
    8) scope: the scope to which a point belongs.
    '''

    def __init__(self,site,rcoord,icoord=[],atom=0,norbital=1,nspin=2,nnambu=1,scope=None):
        self.site=site
        self.rcoord=array(rcoord)
        self.icoord=array(icoord)
        self.atom=atom
        self.norbital=norbital
        self.nspin=nspin
        self.nnambu=nnambu
        self.scope=str(scope)

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        if self.scope=='None':
            return 'Site, atom, norbital, nspin, nnambu: '+str(self.site)+', '+str(self.atom)+', '+str(self.norbital)+', '+str(self.nspin)+', '+str(self.nnambu)+'\n'+'Rcoord: '+str(self.rcoord)+'\n'+'Icoord: '+str(self.icoord)+'\n'
        else:
            return 'Scope, Site, atom, norbital, nspin, nnambu: '+self.scope+', '+str(self.site)+', '+str(self.atom)+', '+str(self.norbital)+', '+str(self.nspin)+', '+str(self.nnambu)+'\n'+'Rcoord: '+str(self.rcoord)+'\n'+'Icoord: '+str(self.icoord)+'\n'

    def __eq__(self,other):
        '''
        Overloaded operator(==).
        '''
        return self.scope==other.scope and self.site==other.site and norm(self.rcoord-other.rcoord)<RZERO and norm(self.icoord-other.icoord)<RZERO and self.atom==other.atom and self.norbital==other.norbital and self.nspin==other.nspin and self.nnambu==other.nnambu
    
    def __ne__(self,other):
        '''
        Overloaded operator(!=).
        '''
        return not self==other
        
    def seq_state(self,orbital,spin,nambu):
        '''
        Sequence of a input state with orbital, spin and nambu index assigned.
        '''
        if nambu in (0,1):
            return spin+orbital*self.nspin+nambu*self.norbital*self.nspin
        else:
            raise ValueError("Point seq_state error: the nambu index must be 0 or 1.")
        
    def state_index(self,seq_state):
        '''
        Orbital, spin and nambu index of a state whose sequence equals the input seq_state.
        '''
        return {'spin':seq_state%(self.norbital*self.nspin)%self.nspin,'orbital':seq_state%(self.norbital*self.nspin)/self.nspin,'nambu':seq_state/(self.norbital*self.nspin)}
