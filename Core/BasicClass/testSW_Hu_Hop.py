# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 16:09:59 2016

@author: DZY
"""
from LatticePy import *
#from DZYHuBbard import *
from DZYSuperbondPy import *

from numpy import kron
def kronecker(*arg):
    n=len(arg)
    out=arg[-1]
    for x in xrange(-2,-n-1,-1):
        out=kron(arg[x],out)
    return out

b=Struct(nspin=2,norbital=3)

p=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=b)

U=5
UP=4
J=1
JP=1

term=HuBbard('U',U,UP,J,JP)

mat=term.mat(p)

out=[]
m=3
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
for i in xrange(m):
    a=zeros((m,m))
    a[i,i]=1
    AP=A-2*kronecker(s[0],a,s[0],a)
BP=A-AP-B
out.append(s3.dot(A-eye(2*m*2*m))/4.0/(UP-JP))
out.append(s1.dot(AP+eye(2*m*2*m))/4.0/(UP+JP))
out.append(s1.dot(BP+B*(m-2.0)/m)/4.0/(U-J))
out.append(s1.dot(2.0/m*B)/4.0/(U+(m-1)*J))

print where(abs(mat[0]-out[0])>1e-10)
print where(abs(mat[1]-out[1])>1e-10)
print where(abs(mat[2]-out[2])>1e-10)
print where(abs(mat[3]-out[3])>1e-10)