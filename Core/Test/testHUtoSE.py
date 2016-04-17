# -*- coding: utf-8 -*-
"""
Created on Mon Apr 04 12:04:49 2016

@author: DZY
"""


from Hamiltonian.Core.BasicClass.LatticePy import *
from Hamiltonian.Core.BasicClass.DZYSuperexchangePy import *
from Hamiltonian.Core.BasicClass.BasisEPy import *
from Hamiltonian.Core.BasicClass.OperatorPy import *
from Hamiltonian.Core.BasicClass.OperatorRepresentationPy import *
from Hamiltonian.Core.BasicAlgorithm.KroneckerPy import *
from Hamiltonian.Core.CoreAlgorithm.ONRPy import *

#from sklearn.preprocessing import normalize
#from scipy.sparse.linalg import eigsh

m=3
b=Fermi(norbital=m,nspin=2)
p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=b)
#p2=Point(site=1,rcoord=[1.0/2.0*sqrt(3),1.0/2.0],icoord=[0.0,0.0],struct=b)
a=Lattice('point',[p1])
c=annihilation(a.table())
basis0=BasisE((6,1))
basis1=BasisE((6,0))
cm0=annihilationRep(c,basis0,basis1)
basis1=BasisE((6,2))
cm2=annihilationRep(c,basis1,basis0)


s=[eye(2, dtype=complex),array([[0,1],[1,0]],dtype=complex),array([[0,-1j],[1j,0]],dtype=complex),array([[1,0],[0,-1]],dtype=complex)]
l=[eye(3, dtype=complex),array([[0,0,0],[0,0,1j],[0,-1j,0]],dtype=complex),array([[0,0,-1j],[0,0,0],[1j,0,0]],dtype=complex),array([[0,1j,0],[-1j,0,0],[0,0,0]],dtype=complex)]

term=SOC(mode='SE',tag='lambda',value=-0.0,bond={'dim':6,'test':0},mat=kron(s[1],l[1])+kron(s[2],l[2])+kron(s[3],l[3]))
Soc=term.to_GUterm(a.table())
#print Soc.mesh(a.bonds[0])
U=5.0
UP=3.0
J=1.0
JP=1.0


Hu=ONR(
        name=       'Hubbard',
        ensemble=   'c',
        filling=    0.0,
        mu=         0.0,
        basis=      BasisE((2*m,2)),
        #basis=      BasisE((2*m*n,m*n)),
        nspin=      1,
        lattice=    a,
        terms=[     Hubbard('U',[U,UP,J]),
                    Soc
                    ]
    )
p1j=Hu.set_projector()


Hu=ONR(
        name=       'Hubbard',
        ensemble=   'c',
        filling=    0.0,
        mu=         0.0,
        basis=      BasisE((2*m,1)),
        #basis=      BasisE((2*m*n,m*n)),
        nspin=      1,
        lattice=    a,
        terms=[     Hubbard('U',[U,UP,J]),
                    Soc
                    ]
    )
p0=Hu.set_projector()

Hu=ONR(
        name=       'Hubbard',
        ensemble=   'c',
        filling=    0.0,
        mu=         0.0,
        basis=      BasisE((2*m,0)),
        #basis=      BasisE((2*m*n,m*n)),
        nspin=      1,
        lattice=    a,
        terms=[     Hubbard('U',[U,UP,J]),
                    Soc
                    ]
    )
p1i=Hu.set_projector()
#
t=kron(array([[1,0],[0,1]]),array([[1,0,0],[0,1,0],[0,0,1]]))
SE=HUtoSE(t,cm0,cm2,p0,p0,p1i,p1j,'SE',tag='J',bond={'dim':36,'test':1})
#
out=[]

l0=eye(m, dtype=complex)
s=[eye(2, dtype=complex),array([[0,1],[1,0]],dtype=complex),array([[0,-1j],[1j,0]],dtype=complex),array([[1,0],[0,-1]],dtype=complex)]
s1=-kronecker(s[0],l0,s[0],l0)+kronecker(s[1],l0,s[1],l0)+kronecker(s[2],l0,s[2],l0)+kronecker(s[3],l0,s[3],l0)
s3=3.0*kronecker(s[0],l0,s[0],l0)+kronecker(s[1],l0,s[1],l0)+kronecker(s[2],l0,s[2],l0)+kronecker(s[3],l0,s[3],l0)
A=zeros((2*m*2*m,2*m*2*m))
AP=zeros((2*m*2*m,2*m*2*m))
B=zeros((2*m*2*m,2*m*2*m))
for i in xrange(m):
    for j in xrange(m):
        a=zeros((m,m))
        a[i,j]=1
        A=A+kronecker(s[0],a,s[0],a.T)
        B=B+kronecker(s[0],a,s[0],a)
AP=A
for i in xrange(m):
    a=zeros((m,m))
    a[i,i]=1
    AP=AP-2*kronecker(s[0],a,s[0],a)
BP=A-AP-B
out.append(s3.dot(A-eye(2*m*2*m))/4.0/(UP-JP))
out.append(s1.dot(AP+eye(2*m*2*m))/4.0/(UP+JP))
out.append(s1.dot(BP+B*(m-2.0)/m)/4.0/(U-J))
out.append(s1.dot(2.0/m*B)/4.0/(U+(m-1)*J))

print where(abs(SE.mat-sum(out,0))>10**-10)