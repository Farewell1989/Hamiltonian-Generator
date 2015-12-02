'''
Basic geometry, including
1) functions: azimuthd,azimuth,polard,polar,volume,is_parallel
2) classes: Point
'''
from StructPy import *
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
        rcoord: 1D ndarray
            The coordinate in real space.
        icoord: 1D ndarray
            the coordinate in lattice space.
        struct: Struct
            The inner structure of the point.
    '''

    def __init__(self,rcoord,icoord=None,struct=None):
        '''
        Constructor.
        Parameters:
            rcoord: 1D array-like
                The coordinate in real space.
            icoord: 1D array-like,optional
                The coordinate in lattice space.
            struct: Struct
                The inner structure of the point.
        '''
        self.rcoord=array(rcoord)
        self.icoord=array([]) if icoord is None else array(icoord)
        self.struct=struct

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        return str(self.struct)+'\n'+'Rcoord: '+str(self.rcoord)+'\n'+'Icoord: '+str(self.icoord)+'\n'

    def __eq__(self,other):
        '''
        Overloaded operator(==).
        '''
        return self.struct==other.struct and norm(self.rcoord-other.rcoord)<RZERO and norm(self.icoord-other.icoord)<RZERO
    
    def __ne__(self,other):
        '''
        Overloaded operator(!=).
        '''
        return not self==other

    @property
    def tag(self):
        '''
        This method returns a string as the tag of a structured point.
        '''
        return self.struct.scope+str(self.struct.site)

def translation(points,vector):
    '''
    This function returns the translated points.
    Parameters:
        points: list of Point
            The original points.
        vector: 1D ndarray
            The translation vector.
    Returns:
        result: list of Point
            The translated points.
    '''
    result=[]
    for p in points:
        result.append(Point(rcoord=p.rcoord+vector,icoord=deepcopy(p.icoord),struct=deepcopy(p.struct)))
    return result

def rotation(points=None,coords=None,angle=0,axis=None,center=None):
    '''
    This function returns the rotated points or coords.
    Parameters:
        points: list of Point
            The original points.
        coords: list of 1D ndarray
            The original coords.
        angle: float
            The rotated angle
        axis: 1D array-like, optional
            The rotation axis. Default the z-axis.
            Not supported yet.
        center: 1D array-like, optional
            The center of the axis. Defualt the origin.
    Returns:
        result: list of Point/1D ndarray
            The rotated points or coords.
    Note: points and coords cannot be both None or not None.
    '''
    if points is None and coords is None:
        raise ValueError('rotation error: both points and coords are None.')
    if points is not None and coords is not None:
        raise ValueError('rotation error: both points and coords are not None.')
    result=[]
    if center is None: center=0
    m11=cos(angle);m21=-sin(angle);m12=-m21;m22=m11
    m=array([[m11,m12],[m21,m22]])
    if points is not None:
        for p in points:
            result.append(Point(rcoord=dot(m,p.rcoord-center)+center,icoord=deepcopy(p.icoord),struct=deepcopy(p.struct)))
    if coords is not None:
        for coord in coords:
            result.append(dot(m,coord-center)+center)
    return result
