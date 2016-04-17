# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 20:14:31 2015

@author: DZY
"""

from LatticePy import *
from DZYSuperexchangePy import *
from DZYSEGeneratorPy import *
from DZYSuperbondPy import *
from DZYSindexPackPy import *

b=Struct(nspin=2)

p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=b)
p2=Point(site=1,rcoord=[1.0,0.0],icoord=[0.0,0.0],struct=b)
p3=Point(site=2,rcoord=[0.0,1.0],icoord=[0.0,0.0],struct=b)
p4=Point(site=3,rcoord=[1.0,1.0],icoord=[0.0,0.0],struct=b)
a=Lattice('link',[p1,p2,p3,p4],nneighbour=2,superbond=2)

pack=(Gmma(2,(0,1))*Gmma(2,(1,0))+Gmma(2,(1,0))*Gmma(2,(0,1)))*2.0+Gmma(2,(0,0))*Gmma(2,(0,0))+Gmma(2,(1,1))*Gmma(2,(1,1))+(Gmma(2,(0,0))*Gmma(2,(1,1))+Gmma(2,(1,1))*Gmma(2,(0,0)))*(-1.0)

term=Superexchange(mode='SE',tag='J',value=1,bond={'dim':4,'test':1},mat=None,sindexpack=pack)
term1=Superexchange(mode='SE',tag='J1',value=0.5,bond={'dim':4,'test':2},mat=None,sindexpack=pack)
term0=Superexchange(mode='SE',tag='J0',value=2,bond={'dim':2,'test':0},mat=array([[1,0],[0,-1]])+array([[0,1],[1,0]])+array([[0,-1j],[1j,0]]),sindexpack=None)

terms=term0+term+term1
H=SEGenerator(a.dim,terms.bondclassify(a.superbonds),terms=terms)



def link(dim,x,y):
    s=[]
    s.append(array([[1,0],[0,1]]))
    s.append(array([[0,1],[1,0]]))
    s.append(array([[0,-1j],[1j,0]]))
    s.append(array([[1,0],[0,-1]]))
    
    out=[1.0,1.0,1.0]
    for i in xrange(dim):
        if i in [x,y]:
            out[0]=sparse.kron(s[1],out[0])
            out[1]=sparse.kron(s[2],out[1])
            out[2]=sparse.kron(s[3],out[2])
        else:
            out[0]=sparse.kron(s[0],out[0])
            out[1]=sparse.kron(s[0],out[1])
            out[2]=sparse.kron(s[0],out[2])
    return out[0]+out[1]+out[2]

hh=link(4,0,1)+link(4,0,2)+link(4,1,3)+link(4,2,3)+0.5*(link(4,0,3)+link(4,1,2))+2*(link(4,1,1)+link(4,0,0)+link(4,2,2)+link(4,3,3))

print (H.hamiltonian-hh).toarray()