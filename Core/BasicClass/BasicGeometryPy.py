'''
Basic geometry, including
1) functions: azimuthd,azimuth,polard,polar,volume,is_parallel
2) classes: Point
'''
from numpy import *
from numpy.linalg import norm
from ConstantPy import RZERO

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
    Structured point.
    Attributes:
        site: integer 
            The site index of the point, start with 0.
        rcoord: 1D ndarray
            The coordinate in real space.
        icoord: 1D ndarray
            the coordinate in lattice space.
        atom: integer, default value 0
            The atom species on this point.
        norbital: integer, default value 1
            Number of orbitals.
        nspin: integer, default value 2
            Number of spins.
        nnambu: integer, default value 1.
            An integer to indicate whether or not using the Nambu space. 1 means no and 2 means yes.
        scope: string, default value 'None'
            The scope to which a point belongs.
    '''

    def __init__(self,site,rcoord,icoord=None,atom=0,norbital=1,nspin=2,nnambu=1,scope=None):
        '''
        Constructor.
        Parameters:
            site: integer 
                The site index of the point.
            rcoord: 1D array-like
                The coordinate in real space.
            icoord: 1D array-like,optional
                the coordinate in lattice space.
            atom: integer, optional
                The atom species on this point.
            norbital: integer, optional
                Number of orbitals.
            nspin: integer, optional
                Number of spins.
            nnambu: integer, optional.
                A flag to indicate whether or not using the Nambu space. 1 means no and 2 means yes.
            scope: string, optional
                The scope to which a point belongs.
        '''
        self.site=site
        self.rcoord=array(rcoord)
        self.icoord=array([]) if icoord is None else array(icoord)
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
        This methods returns the sequence of a input state with orbital, spin and nambu index assigned.
        '''
        if nambu in (0,1):
            return spin+orbital*self.nspin+nambu*self.norbital*self.nspin
        else:
            raise ValueError("Point seq_state error: the nambu index must be 0 or 1.")
        
    def state_index(self,seq_state):
        '''
        This methods returns the orbital, spin and nambu index of a state whose sequence equals the input seq_state.
        Parameters:
            seq_state: integer
                The sequence of the state.
        Returns:
            A dict in the form {'spin':...,'orbital':...,'nambu':...}
        '''
        return {'spin':seq_state%(self.norbital*self.nspin)%self.nspin,'orbital':seq_state%(self.norbital*self.nspin)/self.nspin,'nambu':seq_state/(self.norbital*self.nspin)}

def translation(points,vector,scope=None):
    '''
    This function returns the translated points.
    Parameters:
        points: list of Point
            The original points.
        vector: 1D ndarray
            The translation vector.
        scope: string, optional
            The scope of the translated points.
            When it is None, the translated points share the same scope with the original ones'.
    Returns:
        result: list of Point
            The translated points.
    '''
    result=[]
    for p in points:
        result.append(
            Point(
                scope=scope if scope is None else p.scope,
                site=p.site,
                rcoord=p.rcoord+vector,
                icoord=p.icoord,
                atom=p.atom,
                norbital=p.norbital,
                nspin=p.nspin,
                nnambu=p.nnambu
                )
        )
    return result

def rotation(points,angle,axis=None,center=None,scope=None):
    '''
    This function returns the rotated points.
    Parameters:
        points: list of Point
            The original points.
        angle: float
            The rotated angle
        axis: 1D array-like, optional
            The rotation axis. Default the z-axis.
            Not supported yet.
        center: 1D array-like, optional
            The center of the axis. Defualt the origin.
        scope: string, optional
            The scope of the rotated points.
            When it is None, the rotated points share the same scope with the original ones'.
    Returns:
        result: list of Point
            The rotated points.
    '''
    result=[]
    if center is None: center=0
    for p in points:
        m11=cos(angle);m21=-sin(angle);m12=-m21;m22=m11
        rcoord=dot(array([[m11,m12],[m21,m22]]),p.rcoord-center)+center
        result.append(
            Point(
                scope=scope if scope is None else p.scope,
                site=p.site,
                rcoord=rcoord,
                icoord=p.icoord,
                atom=p.atom,
                norbital=p.norbital,
                nspin=p.nspin,
                nnambu=p.nnambu
                )
        )
    return result
