# -*- coding: utf-8 -*-
"""
Created on Wed Jan 06 19:28:26 2016

@author: DZY
"""


from Hamiltonian.Core.CoreAlgorithm.SpinWPy import *
from Hamiltonian.Core.BasicClass.BaseSpacePy import *
from Hamiltonian.Core.BasicClass.LatticePy import *

J=-1.0
K=10.5/4.0

b=Struct(nspin=2)

p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=b)
p2=Point(site=1,rcoord=[1.0/2.0*sqrt(3),1.0/2.0],icoord=[0.0,0.0],struct=b)
p3=Point(site=2,rcoord=[1.0/2.0*sqrt(3),3.0/2.0],icoord=[0.0,0.0],struct=b)
p4=Point(site=3,rcoord=[1.0*sqrt(3),2.0],icoord=[0.0,0.0],struct=b)
a=Lattice('zigzag',[p1,p2,p3,p4],vectors=[array([1.0*sqrt(3),0.0]),array([1.0*sqrt(3),3.0])],nneighbour=1,superbond=1)

pack=(Gmma(2,(1,0))*Gmma(2,(0,1))+Gmma(2,(0,1))*Gmma(2,(1,0)))*2.0+Gmma(2,(0,0))*Gmma(2,(0,0))+Gmma(2,(1,1))*Gmma(2,(1,1))+(Gmma(2,(0,0))*Gmma(2,(1,1))+Gmma(2,(1,1))*Gmma(2,(0,0)))*(-1.0)

sx=array([[0.0,1.0],[1.0,0.0]])
sy=array([[0.0,-1.0j],[1.0j,0.0]])
sz=array([[1.0,0.0],[0.0,-1.0]])
XX=kron(sx,sx)
YY=kron(sy,sy)
ZZ=kron(sz,sz)
xlink=array([1.0/2.0*sqrt(3),1.0/2.0])
ylink=array([1.0/2.0*sqrt(3),-1.0/2.0])
zlink=array([0.0,1.0])

term=Superexchange(mode='SE',tag='J',value=J,bond={'dim':4,'test':1},mat=None,sindexpack=pack)
termx=Superexchange(mode='SE',tag='Kx',value=2.0*K,bond={'dim':4,'test':lambda bond:(bond.n==2)&is_parallel(bond.rcoords[0]-bond.rcoords[-1],xlink)},mat=XX)
termy=Superexchange(mode='SE',tag='Ky',value=2.0*K,bond={'dim':4,'test':lambda bond:(bond.n==2)&is_parallel(bond.rcoords[0]-bond.rcoords[-1],ylink)},mat=YY)
termz=Superexchange(mode='SE',tag='Kz',value=2.0*K,bond={'dim':4,'test':lambda bond:(bond.n==2)&is_parallel(bond.rcoords[0]-bond.rcoords[-1],zlink)},mat=ZZ)


zigzag=order(a.dim)
zigzag.set_order(1,0,array([[1,0],[0,1]]))
zigzag.set_order(2,0,expm(1j*pi/2*array([[0,-1j],[1j,0]])))
zigzag.set_order(3,0,expm(1j*pi/2*array([[0,-1j],[1j,0]])))

zigzag.settle()

#H=SpinWGenerator(a.dim,zigzag,a.superbonds,terms=term+termx+termy+termz)

Test=SpinW(
    name=       'DZY',
    lattice=    a,
    order=      zigzag,
    terms=      term+termx+termy+termz,
    )
#a.addapps('EB',EB(save_data=False,run=TBAEB))
#a.addapps('DOS',DOS(ne=400,eta=0.01,save_data=False,run=TBADOS))
#a.addapps('EB',EB(path=line_1d(nk=200),save_data=False,run=TBAEB))
#a.addapps('DOS',DOS(BZ=line_1d(nk=10000),eta=0.01,ne=400,save_data=False,run=TBADOS))
Test.addapps('EB',EB(path=rectangle_gym(reciprocals=[array([4.0/3*sqrt(3)*pi,0.0]),array([0.0,4.0/3.0*pi])],nk=100)
,run=TBAEB,save_data=False))
Test.runapps()

