# -*- coding: utf-8 -*-
"""
Created on Mon Feb 29 14:16:10 2016

@author: DZY
"""

from Hamiltonian.Core.CoreAlgorithm.Hubbard_spinWPy import *
from Hamiltonian.Core.BasicClass.BaseSpacePy import *
from Hamiltonian.Core.BasicClass.LatticePy import *

t=1.0
tp=3.0/2
sx=array([[0.0,1.0],[1.0,0.0]])
sy=array([[0.0,-1.0j],[1.0j,0.0]])
sz=array([[1.0,0.0],[0.0,-1.0]])

U=4.0

b=Struct(nspin=2)

p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=b)
p2=Point(site=1,rcoord=[1.0/2.0*sqrt(3),1.0/2.0],icoord=[0.0,0.0],struct=b)
p3=Point(site=2,rcoord=[1.0/2.0*sqrt(3),3.0/2.0],icoord=[0.0,0.0],struct=b)
p4=Point(site=3,rcoord=[1.0*sqrt(3),2.0],icoord=[0.0,0.0],struct=b)
a=Lattice('zigzag',[p1,p2,p3,p4],vectors=[array([1.0*sqrt(3),0.0]),array([1.0*sqrt(3),3.0])],nneighbour=1,superbond=1)



xlink=array([1.0/2.0*sqrt(3),1.0/2.0])
ylink=array([1.0/2.0*sqrt(3),-1.0/2.0])
zlink=array([0.0,1.0])

term=HoPping(mode='HP',tag='t',value=t,bond=lambda bond:bond.tag==1,t=eye(2))
termx=HoPping(mode='HP',tag='tx',value=tp,bond=lambda bond:is_parallel(bond.rcoords[0]-bond.rcoords[-1],xlink) if bond.tag==1 else 0,t=sx)
termy=HoPping(mode='HP',tag='ty',value=tp,bond=lambda bond:is_parallel(bond.rcoords[0]-bond.rcoords[-1],ylink) if bond.tag==1 else 0,t=sy)
termz=HoPping(mode='HP',tag='tz',value=tp,bond=lambda bond:is_parallel(bond.rcoords[0]-bond.rcoords[-1],zlink) if bond.tag==1 else 0,t=sz)

termU=HuBbard(tag='U',U=U)

zigzag=order(a.dim)
zigzag.set_order(1,0,array([[1,0],[0,1]]))
zigzag.set_order(2,0,expm(1j*pi/2*array([[0,-1j],[1j,0]])))
zigzag.set_order(3,0,expm(1j*pi/2*array([[0,-1j],[1j,0]])))

zigzag.settle()

#H=Hubbard_SpinWGenerator(a.dim,zigzag,a.superbonds,hopping=[term,termx,termy,termz],hubbard=[termU])

Test=Hubbard_SpinW(
    name=       'DZY',
    lattice=    a,
    order=      zigzag,
    hopping=    [term,termx,termy,termz],
    hubbard=    [termU],
    terms=      [],
    )
#a.addapps('EB',EB(save_data=False,run=TBAEB))
#a.addapps('DOS',DOS(ne=400,eta=0.01,save_data=False,run=TBADOS))
#a.addapps('EB',EB(path=line_1d(nk=200),save_data=False,run=TBAEB))
#a.addapps('DOS',DOS(BZ=line_1d(nk=10000),eta=0.01,ne=400,save_data=False,run=TBADOS))
Test.addapps('EB',EB(path=rectangle_gym(reciprocals=[array([4.0/3*sqrt(3)*pi,0.0]),array([0.0,4.0/3.0*pi])],nk=100)
,run=TBAEB,save_data=False))
Test.runapps()