from BasicGeometryPy import *
from numpy.linalg import inv
import matplotlib.pyplot as plt
import time
class BaseSpace:
    '''
    The BaseSpace class provides a unified description of all kinds of parameter spaces. It has the following attributes:
    1) mode: a string specifies the meaning of the coordinate in the base space.
    2) mesh: a two-dimensional numpy array, of whose shape, the first element equals the number of points in the base space and the second element equals the dimension of the base space.
    3) volume: the volume of the base space.
    '''
    def __init__(self,mode,mesh=None,volume=0.0):
        self.mode=mode
        self.mesh=mesh
        self.volume=volume

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        return str(self.mesh)

    @property
    def rank(self):
        '''
        The total number of points in the base space, i.e. self.mesh.shape[0].
        '''
        return self.mesh.shape[0]

    def plot(self,show=True,name='BaseSpace'):
        '''
        Plot all the points contained in its mesh. Only two dimensional base spaces are supported.
        '''
        plt.axis('equal')
        x=self.mesh[:,0]
        y=self.mesh[:,1]
        plt.scatter(x,y)
        if show:
            plt.show()
        else:
            plt.savefig(name+'.png')
        plt.close()

def KSpace(reciprocals=None,nk=100,mesh=None,volume=0.0):
    '''
    The KSpace function provides a unified description of the whole Broullouin zone(BZ), a path in the BZ, or just some isolated points in the BZ. There are two ways to initialize a non-empty instance:
    1) Assign the reciprocals and the number of slices nk. The mesh and its volume in K-space is generated automatically.
    2) Assign the mesh directly. In this case, the attribute 'volume' may have no physical meanings.
    '''
    result=BaseSpace(mode='k',mesh=mesh,volume=volume)
    if not reciprocals is None:
        nvectors=len(reciprocals)
        if nvectors==1:
            result.volume=norm(reciprocals[0])
        elif nvectors==2:
            result.volume=abs(cross(reciprocals[0],reciprocals[1]))
        elif nvectors==3:
            result.volume==abs(volume(reciprocals[0],reciprocals[1],reciprocals[2]))
        else:
            raise ValueError("KSpace error: the number of reciprocals should not be greater than 3.")
        ndim=reciprocals[0].shape[0]
        ubi=1;ubj=1;ubk=1
        if nvectors>=1:ubi=nk
        if nvectors>=2:ubj=nk
        if nvectors>=3:ubk=nk
        result.mesh=zeros((ubi*ubj*ubk,ndim))
        for i in xrange(ubi):
            for j in xrange(ubj):
                for k in xrange(ubk):
                    for l in xrange(ndim):
                        for h in xrange(nvectors):
                            if h==0: buff=1.0*i/ubi-0.5
                            if h==1: buff=1.0*j/ubj-0.5
                            if h==2: buff=1.0*k/ubk-0.5
                            result.mesh[(i*ubj+j)*ubk+k,l]=result.mesh[(i*ubj+j)*ubk+k,l]+reciprocals[h][l]*buff
    return result

def line_1d(reciprocals=None,nk=100):
    '''
    The BZ of 1D K-space.
    '''
    if reciprocals is None:
        recips=[array([2*pi])]
    else:
        recips=reciprocals
    return KSpace(reciprocals=recips,nk=nk)

def rectangle_gxm(reciprocals=None,nk=100):
    '''
    The Gamma-X-M-Gamma path in the rectangular BZ.
    '''
    result=KSpace(nk=nk)
    if not reciprocals is None:
        b1=reciprocals[0]
        b2=reciprocals[1]
    else:
        b1=array([2*pi,0.0])
        b2=array([0.0,2*pi])
    ndim=b1.shape[0]
    result.mesh=zeros((3*nk,ndim))
    for i in xrange(nk):
        result.mesh[i,:]=b1/2*i/nk
        result.mesh[nk+i,:]=b1/2+b2/2*i/nk
        result.mesh[nk*2+i,:]=(b1+b2)/2*(1-1.0*i/nk)
    return result

def rectangle_gym(reciprocals=None,nk=100):
    '''
    The Gamma-X-M-Gamma path in the rectangular BZ.
    '''
    result=KSpace(nk=nk)
    if not reciprocals is None:
        b1=reciprocals[1]
        b2=reciprocals[0]
    else:
        b1=array([0.0,2*pi])
        b2=array([2*pi,0.0])
    ndim=b1.shape[0]
    result.mesh=zeros((3*nk,ndim))
    for i in xrange(nk):
        result.mesh[i,:]=b1/2*i/self.nk
        result.mesh[nk+i,:]=b1/2+b2/2*i/nk
        result.mesh[nk*2+i,:]=(b1+b2)/2*(1-1.0*i/nk)
    return result

def rectangle_bz(reciprocals=None,nk=100):
    '''
    The whole rectangular BZ.
    '''
    if reciprocals is None:
        recips=[]
        recips.append(array([2*pi,0.0]))
        recips.append(array([0.0,2*pi]))
    else:
        recips=reciprocals
    return KSpace(reciprocals=recips,nk=nk)

square_gxm=rectangle_gxm
square_bz=rectangle_bz

def hexagon_gkm(reciprocals=None,nk=100,vh='H'):
    '''
    The Gamma-K-M-Gamma path in the hexagonal BZ.
    '''
    result=KSpace(nk=nk)
    if not reciprocals is None:
        b1=reciprocals[0]
        b2=reciprocals[1]
        if abs(inner(b1,b2)+0.5)<RZERO:
            b2=-b2
        elif abs(inner(b1,b2)-0.5)>RZERO:
            raise ValueError("Hexagon_gkm error: the reciprocals are too wired.")
    else:
        if vh in ('H','h'):
          b1=array([sqrt(3.0)/2,0.5])*4*pi/sqrt(3.0)
          b2=array([sqrt(3.0)/2,-0.5])*4*pi/sqrt(3.0)
        else:
          b1=array([1.0,0.0])*4*pi/sqrt(3.0)
          b2=array([0.5,sqrt(3.0)/2])*4*pi/sqrt(3.0)
    ndim=b1.shape[0]
    result.mesh=zeros((3*nk,ndim))
    for i in xrange(nk):
        result.mesh[i,:]=(b1+b2)/3*i/nk
        result.mesh[nk+i,:]=(b1+b2)/3+(b1-2*b2)/6*i/nk
        result.mesh[nk*2+i,:]=b1/2*(1-1.0*i/nk)
    return result

def hexagon_bz(reciprocals=None,nk=100,vh='H'):
    '''
    The whole hexagonal BZ.
    '''
    result=KSpace(nk=nk)
    if not reciprocals is None:
        b1=reciprocals[0]
        b2=reciprocals[1]
        if abs(inner(b1,b2)+0.5)<RZERO:
            b2=-b2
        elif abs(inner(b1,b2)-0.5)>RZERO:
            raise ValueError("Hexagon_gkm error: the reciprocals are too wired.")
    else:
        if vh in ('H','h'):
          b1=array([sqrt(3.0)/2,0.5])*4*pi/sqrt(3.0)
          b2=array([sqrt(3.0)/2,-0.5])*4*pi/sqrt(3.0)
        else:
          b1=array([1.0,0.0])*4*pi/sqrt(3.0)
          b2=array([0.5,sqrt(3.0)/2])*4*pi/sqrt(3.0)
    ndim=b1.shape[0]
    result.mesh=zeros((nk**2,ndim))
    p0=-(b1+b2)/3
    p1=(b1+b2)/3
    p2=(b1+b2)*2/3
    p3=(b1*2-b2)/3
    p4=(b2*2-b1)/3
    for i in xrange(nk):
        for j in xrange(nk):
          coords=b1*(i-1)/nk+b2*(j-1)/nk+p0
          if in_triangle(coords,p1,p2,p3): coords=coords-b1
          if in_triangle(coords,p1,p2,p4): coords=coords-b2
          result.mesh[i*nk+j,:]=coords
    result.volume=abs(cross(b1,b2))
    return result
    
def in_triangle(p0,p1,p2,p3):
    '''
    Judge whether a point represented by p0 belongs to the interior of a triangle whose vertices are p1,p2 and p3.
    '''
    a=zeros((3,3))
    b=zeros(3)
    x=zeros(3)
    ndim=p0.shape[0]
    a[0:ndim,0]=p2-p1
    a[0:ndim,1]=p3-p1
    a[(2 if ndim==2 else 0):3,2]=cross(p2-p1,p3-p2)
    b[0:ndim]=p0-p1
    x=dot(inv(a),b)
    if x[0]>=0 and x[0]<=1 and x[1]>=0 and x[1]<=1 and x[0]+x[1]<=1:
        return True
    else:
        return False

def TSpace(mesh):
    '''
    The time space.
    '''
    return BaseSpace(mode='t',mesh=mesh,volume=mesh.max()-mesh.min())

# The following codes are used for tests only.
def test_basespace():
    a=KSpace(reciprocals=[array([2*pi,0.0]),array([0.0,2*pi])],nk=100)
    a.plot(show=True)
    print a.volume/(2*pi)**2
    a=KSpace(reciprocals=[array([1.0,0.0]),array([0.5,sqrt(3.0)/2])],nk=100)
    a.plot(show=True)
    print a.volume
    square_gxm(nk=100).plot()

def test_basespace_functions():
    a=square_bz(reciprocals=[array([1.0,1.0]),array([1.0,-1.0])],nk=100)
    a.plot(show=True)
    print a.volume
    a=rectangle_bz(nk=100)
    a.plot(show=True)
    print a.volume/(2*pi)**2
    a=hexagon_bz(nk=100,vh='v')
    a.plot(show=True)
    print a.volume/(2*pi)**2