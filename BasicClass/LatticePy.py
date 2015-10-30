from IndexPy import *
from BasicGeometryPy import *
from BondPy import *
from numpy.linalg import inv
from itertools import product
import matplotlib.pyplot as plt

class Lattice:
    '''
    The lattice class provides a unified description of 1D, quasi 1D, 2D, quasi 2D and 3D lattice systems. It has the following attributes:
    1) name: a string representing the lattice's name;
    2) points: a list representing the lattice points in a unit cell;
    3) vectors: a list representing the translation vectors;
    4) reciprocals: a list representing the dual translation vectors;
    5) nneighbour: a number representing the highest order of neighbours;
    6) bonds: a nested list whose first element represent the zero-th order bonds(i.e. the intra unit cell points), second element the nearest neighbour bonds, third element the next-nearest neighbour bonds, etc.
    7) priority: a string to indicate the sequence priority of the allowed indices that can be defined on this lattice with the default value 'NSCO', where 'N','S','C','O' stands for 'nambu', 'spin', 'site' and 'orbital' respectively.
    '''
    
    def __init__(self,name,points,translations=None,vectors=[],nneighbour=1,priority='NSCO'):
        '''
        There are 2 ways to initialize a lattice instance:
        1) Assign the points in a unit cell and the translation vectors.
        2) Assign the points in a small cluster and the translation tuples in the form ((a1,m1),(a2,m2),...), where ai denotes a translation vector and mi is the number of slices along that direction. The corresponding super unit cell will be created and its translation vectors are assigned by the parameter vectors.
        In both cases, the nneighbour attribute can be explicitly assigned while the reciprocals and bonds attributes will be automatically generated.
        '''
        self.name=name
        self.points=points
        if not translations is None:
            for (a,m) in translations:
                npoints=len(self.points)
                for point in self.points[:]:
                    for i in xrange(1,m):
                        self.points.append(Point(site=point.site+i*npoints,rcoord=point.rcoord+i*a,icoord=point.icoord,atom=point.atom,norbital=point.norbital,nspin=point.nspin,nnambu=point.nnambu,scope=point.scope))
            self.points.sort(key=lambda point : point.site)
        self.vectors=vectors
        self.reciprocals=reciprocals(self.vectors)
        self.nneighbour=nneighbour
        self.bonds=bonds(self.points,self.vectors,self.nneighbour)
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
        for nb,bonds in enumerate(self.bonds):
            if nb==1: color='k'
            elif nb==2: color='r'
            elif nb==3: color='b'
            else: color=str(nb*1.0/self.nneighbour)
            for bond in bonds:
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

    def indices(self,nambu=False):
        '''
        Return all the allowed indices that can be defined on this lattice.
        '''
        result=[]
        for point in self.points:
            for orbital in xrange(point.norbital):
                for spin in xrange(point.nspin):
                    if nambu:
                        for buff in xrange(point.nnambu):
                            result.append(Index(scope=self.name,site=point.site,orbital=orbital,spin=spin,nambu=buff))
                    else:
                        result.append(Index(scope=self.name,site=point.site,orbital=orbital,spin=spin,nambu=ANNIHILATION))
        buff=[]
        for i in result:
            buff.append(i.to_str(indication='p'+self.priority))
        result=[]
        for i in sorted(buff):
            result.append(to_index(str=i,indication='p'+self.priority))
        return result

def bonds(points,vectors=None,nneighbour=1):
    '''
    Return all the bonds up to the nneighbour-th neighbour within the input cluster whose points are assigned. If the vectors are non-empty, the inter cluster bonds will be included too.
    Note: the input points will be included in the returned list as the zero-th neighbour bonds.
    '''
    nvectors=len(vectors)
    npoints=len(points)
    ndim=points[0].rcoord.shape[0]
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
    indices=list(product(xrange(-sup1,sup1+1),xrange(-sup2,sup2+1),xrange(-sup3,sup3+1)))
    for i,j,k in indices:
        indices.remove((-i,-j,-k))
        if nvectors==0: disp=zeros(ndim)
        if nvectors>=1: disp=vectors[0]*i
        if nvectors>=2: disp+=vectors[1]*j
        if nvectors>=3: disp+=vectors[2]*k
        for m in xrange(npoints):
            for n in xrange(npoints):
                if i==0 and j==0 and k==0 and n>m: continue
                coord1=points[m].rcoord
                coord2=points[n].rcoord+disp
                cdist=norm(coord1-coord2)
                for l,dist in enumerate(mdists[:]):
                    if abs(cdist)<RZERO:
                        result[0].append(Bond(0,points[m],points[n]))
                        break
                    elif abs(cdist-dist)<RZERO:
                        result[l+1].append(Bond(l+1,points[m],Point(site=points[n].site,rcoord=coord2,icoord=disp,atom=points[n].atom,norbital=points[n].norbital,nspin=points[n].nspin,nnambu=points[n].nnambu,scope=points[n].scope)))
                        break
                    elif cdist<dist:
                        mdists.insert(l,cdist)
                        result.insert(l+1,[])
                        result[l+1].append(Bond(l+1,points[m],Point(site=points[n].site,rcoord=coord2,icoord=disp,atom=points[n].atom,norbital=points[n].norbital,nspin=points[n].nspin,nnambu=points[n].nnambu,scope=points[n].scope)))
                        break
                else:
                    if cdist<RZERO:
                        result[0].append(Bond(0,points[m],points[n]))
                    elif len(mdists)<nneighbour:
                        mdists.append(cdist)
                        result.append([])
                        result[len(result)-1].append(Bond(len(result),points[m],Point(site=points[n].site,rcoord=coord2,icoord=disp,atom=points[n].atom,norbital=points[n].norbital,nspin=points[n].nspin,nnambu=points[n].nnambu,scope=points[n].scope)))
                if len(mdists)>nneighbour:
                    del mdists[nneighbour]
                    del result[nneighbour+1]
    for nb,bonds in enumerate(result):
        for bond in bonds:
            bond.neighbour=nb
    return result

def reciprocals(vectors):
    '''
    Return the corresponding reciprocals dual to the input vectors.
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
