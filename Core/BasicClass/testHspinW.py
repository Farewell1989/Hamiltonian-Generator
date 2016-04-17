# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 14:28:47 2015

@author: DZY
"""

from numpy import *
from DZYOrderPy import *
from DZYSuperexchangePy import *
from DZYSpinWGeneratorPy import *
from DZYSpinWmatrixPy import *
from DZYSuperbondPy import *
from DZYSindexPackPy import *
from BaseSpacePy import *




# -*- coding: utf-8 -*-
"""
Created on Wed Jan 06 19:28:26 2016

@author: DZY
"""

from Hamiltonian.Core.BasicClass.BaseSpacePy import *
from Hamiltonian.Core.BasicClass.LatticePy import *
from Hamiltonian.Core.CoreAlgorithm.Hubbard_spinWPy import *

#from Hubbard_SpinWGenerator import *

tx=zeros((3,3))
tx[1,1]=1.0
tx[2,2]=1.0
ty=zeros((3,3))
ty[0,0]=1.0
ty[2,2]=1.0
tz=zeros((3,3))
tz[0,0]=1.0
tz[1,1]=1.0

U=1.0
J=2.0
Jp=-6.0
Up=U-J-Jp

b=Struct(nspin=2,norbital=3)

p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=b)
p2=Point(site=1,rcoord=[1.0/2.0*sqrt(3),1.0/2.0],icoord=[0.0,0.0],struct=b)
p3=Point(site=2,rcoord=[1.0/2.0*sqrt(3),3.0/2.0],icoord=[0.0,0.0],struct=b)
p4=Point(site=3,rcoord=[1.0*sqrt(3),2.0],icoord=[0.0,0.0],struct=b)
a=Lattice('zigzag',[p1,p2,p3,p4],vectors=[array([1.0*sqrt(3),0.0]),array([1.0*sqrt(3),3.0])],nneighbour=1,superbond=1)

xlink=array([1.0/2.0*sqrt(3),1.0/2.0])
ylink=array([1.0/2.0*sqrt(3),-1.0/2.0])
zlink=array([0.0,1.0])

termx=HoPping(mode='HP',tag='tx',value=1.0,bond=lambda bond:is_parallel(bond.rcoords[0]-bond.rcoords[-1],xlink) if bond.tag==1 else 0,tspin=eye(2),torbital=tx)
termy=HoPping(mode='HP',tag='ty',value=1.0,bond=lambda bond:is_parallel(bond.rcoords[0]-bond.rcoords[-1],ylink) if bond.tag==1 else 0,tspin=eye(2),torbital=ty)
termz=HoPping(mode='HP',tag='tz',value=1.0,bond=lambda bond:is_parallel(bond.rcoords[0]-bond.rcoords[-1],zlink) if bond.tag==1 else 0,tspin=eye(2),torbital=tz)

termU=HuBbard(tag='U',U=U,UP=Up,J=J,JP=Jp)

s=[eye(2, dtype=complex),array([[0,1],[1,0]],dtype=complex),array([[0,-1j],[1j,0]],dtype=complex),array([[1,0],[0,-1]],dtype=complex)]
l=[eye(3, dtype=complex),array([[0,0,0],[0,0,1j],[0,-1j,0]],dtype=complex),array([[0,0,-1j],[0,0,0],[1j,0,0]],dtype=complex),array([[0,1j,0],[-1j,0,0],[0,0,0]],dtype=complex)]

term=Superexchange(mode='SE',tag='lambda',value=-5.0,bond={'dim':6,'test':0},mat=kron(s[1],l[1])+kron(s[2],l[2])+kron(s[3],l[3]))

zigzag=order(a.dim)
zigzag.set_order(1,0,kron(s[0],l[0]))
zigzag.set_order(2,0,expm(1j*pi/2*(kron(s[2]/2.0,l[0])-kron(s[0],l[2]))))
zigzag.set_order(3,0,expm(1j*pi/2*(kron(s[2]/2.0,l[0])-kron(s[0],l[2]))))

zigzag.settle()

#H=Hubbard_SpinWGenerator(a.dim,zigzag,a.superbonds,hopping=[termx,termy,termz],hubbard=[termU],terms=[term])
#
Test=Hubbard_SpinW(
    dout=       'D:\python\picture',
    name=       'DZY',
    lattice=    a,
    order=      zigzag,
    hopping=    [termx,termy,termz],
    hubbard=    [termU],
    terms=      [term],
    )
#a.addapps('EB',EB(save_data=False,run=TBAEB))
#a.addapps('DOS',DOS(ne=400,eta=0.01,save_data=False,run=TBADOS))
#a.addapps('EB',EB(path=line_1d(nk=200),save_data=False,run=TBAEB))
#a.addapps('DOS',DOS(BZ=line_1d(nk=10000),eta=0.01,ne=400,save_data=False,run=TBADOS))
Test.addapps('EB',EB(path=rectangle_gym(reciprocals=[array([4.0/3*sqrt(3)*pi,0.0]),array([0.0,4.0/3.0*pi])],nk=100)
,run=TBAEB,save_data=True))
Test.runapps()


#soc=-kron(s[1],l[1])-kron(s[2],l[2])-kron(s[3],l[3])
#
#sx=kron(s[1]/2.0,l[0])-kron(s[0],l[1])
#sy=kron(s[2]/2.0,l[0])-kron(s[0],l[2])
#sz=kron(s[3]/2.0,l[0])-kron(s[0],l[3])
#
#w,v=linalg.eig(soc+0.1*sz)
#
#V=kron(v[:,(4,0)],v[:,(4,0)])    
#
#test=V.conj().T.dot(H.origin_cache[4].mat+H.origin_cache[4].mat.conj().T).dot(V)
#
#T=zeros((4,4),dtype=complex)
#for x in xrange(4):
#    for y in xrange(4):
#        T[x,y]=4.0*trace(test.dot(kron(s[x],s[y])))
#        
#print where(abs(T)>1e-9),T[where(abs(T)>1e-9)],H.origin_cache[4].bond.sites