'''
Lattice.
'''
from BasicGeometryPy import *
from BondPy import *
from numpy.linalg import inv
import itertools
import matplotlib.pyplot as plt

class Lattice(object):
    '''
    This class provides a unified description of 1D, quasi 1D, 2D, quasi 2D and 3D lattice systems.
    Attributes:
        name: string
            The lattice's name.
        points: dict of Point
            The lattice points in a unit cell.
        vectors: list of 1D ndarray
            The translation vectors.
        reciprocals: list of 1D ndarray
            The dual translation vectors.
        nneighbour: integer
            The highest order of neighbours;
        bonds: list of Bond
            The bonds of the lattice system.
        priority: string
            The sequence priority of the allowed indices that can be defined on this lattice.
            Default value 'PNSCO', where 'P','N','S','C','O' stands for 'scope', 'nambu', 'spin', 'site' and 'orbital' respectively.
    '''
    
    def __init__(self,name,points,translations=None,vectors=[],nneighbour=1,priority='PNSCO'):
        '''
        Constructor.
        It can be used in the following ways:
        1) Lattice(name=...,points=...,translations=...,vectors=...,nneighbour=...,priority=...)
        2) Lattice(name=...,points=...,vectors=...,nneighbour=...,priority=...)
        Parameters:
            name: string
                The name of the lattice.
            points: list of Point
                Two cases:
                1) translations is None: the lattice points in a unit cell.
                2) translations is not None: the original translated points.
            translations: list of 2-tuple, optional
                For each tuple:
                    tuple[0]: 1D ndarray
                        The translation vector for the original points.
                    tuple[1]: integer
                        The number of slices along the corresponding vector.
            vectors: list of 1D ndarray, optional
                The translation vectors of the lattice.
            nneighbour: integer, optional
                The highest order of neighbours.
            priority: string, optional
                The sequence priority of the allowed indices that can be defined on this lattice.
        '''
        self.name=name
        self.points={}
        for p in points:
            if name==p.scope:
                self.points[name+str(p.site)]=p
            else:
                self.points[name+str(p.site)]=Point(scope=name,site=p.site,rcoord=p.rcoord,icoord=p.icoord,struct=p.struct)
        if translations is not None:
            for (a,m) in translations:
                ps=self.points.values()
                inc=max([p.site for p in ps])+1
                for p in ps[:]:
                    for i in xrange(1,m):
                        site=p.site+i*inc
                        self.points[name+str(site)]=Point(scope=p.scope,site=site,rcoord=p.rcoord+i*a,icoord=p.icoord,struct=p.struct)
        self.vectors=vectors
        self.reciprocals=reciprocals(self.vectors)
        self.nneighbour=nneighbour
        self.bonds=[b for bs in bonds(self.points,self.vectors,self.nneighbour) for b in bs]
        self.priority=priority

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result=''
        for bonds in self.bonds:
            for bond in bonds:
                result+=str(bond)
        return result

    def plot(self,show=True):
        '''
        Plot the lattice points and bonds. Only 2D or quasi 1D systems are supported.
        '''
        plt.axes(frameon=0)
        plt.axis('equal')
        plt.title(self.name)
        for bond in self.bonds:
            nb=bond.neighbour
            if nb==1: color='k'
            elif nb==2: color='r'
            elif nb==3: color='b'
            else: color=str(nb*1.0/self.nneighbour)
            if nb==0:
                plt.scatter(bond.spoint.rcoord[0],bond.spoint.rcoord[1])
            else:
                if bond.is_intra_cell():
                    plt.plot([bond.spoint.rcoord[0],bond.epoint.rcoord[0]],[bond.spoint.rcoord[1],bond.epoint.rcoord[1]],color=color)
                else:
                    plt.plot([bond.spoint.rcoord[0],bond.epoint.rcoord[0]],[bond.spoint.rcoord[1],bond.epoint.rcoord[1]],color=color,ls='--')
        frame=plt.gca()
        frame.axes.get_xaxis().set_visible(False)
        frame.axes.get_yaxis().set_visible(False)
        if show:
            plt.show()
        else:
            plt.savefig(self.name+'.png')
        plt.close()

    def table(self,nambu=False):
        '''
        Return a Table instance that contains all the allowed indices which can be defined on this lattice.
        '''
        return union([p.table(nambu=nambu) for p in self.points.itervalues()],key=lambda value: value.to_tuple(indication=self.priority))

def bonds(points,vectors=None,nneighbour=1):
    '''
    This function returns all the bonds up to the nneighbour-th order.
    Parameters:
        points: dict of Point
            The cluster within which the bonds are looked for. 
            If the parameter vectors is not None, the inter cluster bonds will also be searched.
        vectors: list of 1D ndarray, optional
            The translation vectors for the cluster.
        nneighbour: integer, optional
            The highest order of neighbour to be searched.
    Returns:
        result: list of list of Bond
            result[k] contains all the k-th neighbour bonds.
            Note that the input points will be taken as both the start points and end points of the zero-th neighbour bonds.
    '''
    nvectors=len(vectors)
    ndim=points.values()[0].rcoord.shape[0]
    result=[[]]
    if nvectors==0:
        sup1=0;sup2=0;sup3=0
    elif nvectors==1:
        sup1=nneighbour;sup2=0;sup3=0
    elif nvectors==2:
        sup1=nneighbour;sup2=nneighbour;sup3=0
    elif nvectors==3:
        sup1=nneighbour;sup2=nneighbour;sup3=nneighbour
    else:
        raise ValueError("Bonds error: the number of vectors should not be greater than 3.")
    mdists=[]
    indices=list(itertools.product(xrange(-sup1,sup1+1),xrange(-sup2,sup2+1),xrange(-sup3,sup3+1)))
    for i,j,k in indices:
        indices.remove((-i,-j,-k))
        if nvectors==0: disp=zeros(ndim)
        if nvectors>=1: disp=vectors[0]*i
        if nvectors>=2: disp+=vectors[1]*j
        if nvectors>=3: disp+=vectors[2]*k
        for m in points.iterkeys():
            for n in points.iterkeys():
                if i==0 and j==0 and k==0 and n>m: continue
                coord1=points[m].rcoord
                coord2=points[n].rcoord+disp
                cdist=norm(coord1-coord2)
                for l,dist in enumerate(mdists[:]):
                    if abs(cdist)<RZERO:
                        result[0].append(Bond(0,points[m],points[n]))
                        break
                    elif abs(cdist-dist)<RZERO:
                        result[l+1].append(Bond(l+1,points[m],Point(scope=points[n].scope,site=points[n].site,rcoord=coord2,icoord=disp,struct=points[n].struct)))
                        break
                    elif cdist<dist:
                        mdists.insert(l,cdist)
                        result.insert(l+1,[])
                        result[l+1].append(Bond(l+1,points[m],Point(scope=points[n].scope,site=points[n].site,rcoord=coord2,icoord=disp,struct=points[n].struct)))
                        break
                else:
                    if cdist<RZERO:
                        result[0].append(Bond(0,points[m],points[n]))
                    elif len(mdists)<nneighbour:
                        mdists.append(cdist)
                        result.append([])
                        result[len(result)-1].append(Bond(len(result),points[m],Point(scope=points[n].scope,site=points[n].site,rcoord=coord2,icoord=disp,struct=points[n].struct)))
                if len(mdists)>nneighbour:
                    del mdists[nneighbour]
                    del result[nneighbour+1]
    for nb,bonds in enumerate(result):
        for bond in bonds:
            bond.neighbour=nb
    return result

def reciprocals(vectors):
    '''
    This function returns the corresponding reciprocals dual to the input vectors.
    '''
    result=[]
    nvectors=len(vectors)
    if nvectors==0:
        return
    if nvectors==1:
        result.append(array(vectors[0]/(norm(vectors[0]))**2*2*pi))
    elif nvectors in (2,3):
        ndim=vectors[0].shape[0]
        buff=zeros((3,3))
        buff[0:ndim,0]=vectors[0]
        buff[0:ndim,1]=vectors[1]
        if nvectors==2:
            buff[(2 if ndim==2 else 0):3,2]=cross(vectors[0],vectors[1])
        else:
            buff[0:ndim,2]=vectors[2]
        buff=inv(buff)
        result.append(array(buff[0,0:ndim]*2*pi))
        result.append(array(buff[1,0:ndim]*2*pi))
        if nvectors==3:
            result.append(array(buff[2,0:ndim]*2*pi))
    else:
        raise ValueError('Reciprocals error: the number of translation vectors should not be greater than 3.')
    return result

def SuperLattice(name,sublattices,vectors=[],nneighbour=1,priority='PNSCO'):
    '''
    This function returns the union of sublattices.
    Parameters:
        name: string
            The name of the super-lattice.
        sublattices: list of Lattice
            The sub-lattices of the super-lattice.
        vectors: list of 1D ndarray, optional
            The translation vectors of the super-lattice.
        nneighbour: integer,optional
            The highest order of neighbours.
        priority: string, optional
            The sequence priority of the allowed indices that can be defined on this super-lattice.
    Returns:
        result: Lattice
            The super-lattice.
    '''
    result=object.__new__(Lattice)
    result.name=name
    result.points={}
    for lattice in sublattices:
        result.points.update(lattice.points)
    result.vectors=vectors
    result.reciprocals=reciprocals(vectors)
    result.nneighbour=nneighbour
    result.bonds=[b for bs in bonds(result.points,vectors,nneighbour) for b in bs]
    result.priority=priority
    result.sublattices=sublattices
    return result
