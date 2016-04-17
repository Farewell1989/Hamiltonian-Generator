# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 20:09:07 2016

@author: DZY
"""

'''
HuBbard (onsite) interaction terms.
'''

from TermPy import *
from DZYSuperbondPy import *
from BasicGeometryPy import *
from DZYSpinWmatrixPy import *
from Hamiltonian.Core.BasicAlgorithm.KroneckerPy import *

class HuBbard(Term):
    '''
    '''

    def __init__(self,tag,U,UP=0,J=0,JP=0,point=None,modulate=None):
        '''
        Constructor.
        '''
        super(HuBbard,self).__init__('hb',tag,{'U':U,'UP':UP,'J':J,'JP':JP},modulate)
        if point is not None: self.point=point

#    def __str__(self):
#        '''
#        Convert an instance to string.
#        '''
#        result=''
#        if hasattr(self,'atom'): result+='Atom: '+str(self.atom)+'\n'
#        result+='Tag,value: '+self.tag+','+str(self.value)+'\n'
#        return result

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a HuBbard instance with a HuBbard/HuBbardList instance.
        '''
        result=HuBbardList(deepcopy(self))
        if isinstance(other,HuBbard):
            result.append(deepcopy(other))
        elif isinstance(other,HuBbardList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('HuBbard "+" error: the other parameter must be an instance of HuBbard or HuBbardList')
        return result

    def __pos__(self):
        '''
        Overloaded operator(+), i.e. +self.
        '''
        result=HuBbardList(deepcopy(self))
        return result
    
    def mat(self,point):
        out=[]
        if point.struct.nspin!=2:
            raise ValueError('HuBbard error: The point.nspin must be 2')
        s=[eye(2, dtype=complex),array([[0,1],[1,0]],dtype=complex),array([[0,-1j],[1j,0]],dtype=complex),array([[1,0],[0,-1]],dtype=complex)]
        if hasattr(point.struct,'norbital'):
            m=point.struct.norbital
            l0=eye(m, dtype=complex)
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
            out.append(s3.dot(A-eye(2*m*2*m))/4.0/(self.value['UP']-self.value['JP']))
            out.append(s1.dot(AP+eye(2*m*2*m))/4.0/(self.value['UP']+self.value['JP']))
            out.append(s1.dot(BP+B*(m-2.0)/m)/4.0/(self.value['U']-self.value['J']))
            out.append(s1.dot(2.0/m*B)/4.0/(self.value['U']+(m-1)*self.value['J']))
        else:
            out.append((-kron(s[0],s[0])+kron(s[1],s[1])+kron(s[2],s[2])+kron(s[3],s[3]))/2.0/self.value['U'])
        return out
        
    def spinWmat(self,bond):
        out=SpinWmatrix(bond,sum(self.mat(bond.points[0]),0))
        result=triu(out.mat,1)
        result=result+diag(diag(out.mat)/2.0)
        out.mat=result
        out.ormat=out.mat
        return out


class HuBbardList(list):
    '''
    This class pack several HuBbard instances as a whole for convenience.
    '''
    def __init__(self,*arg):
        self.mode='hb'
        for obj in arg:
            if isinstance(obj,HuBbard):
                self.append(obj)
            else:
                raise ValueError("HuBbardList init error: the input argument should be HuBbard instances.")

#    def __str__(self):
#        '''
#        Convert an instance to string.
#        '''
#        result='HuBbard terms:\n'
#        for i,v in enumerate(self):
#            result+='Node '+str(i)+':\n'+str(v)

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a HuBbardList instance with a HuBbard/HuBbardList instance.
        '''
        result=HuBbardList(*deepcopy(self))
        if isinstance(other,HuBbard):
            result.append(deepcopy(other))
        elif isinstance(other,HuBbardList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('HuBbardList "+" error: the other parameter must be an instance of HuBbard or HuBbardList')
        return result

    def __mul__(self,other):
        '''
        Overloaded operator(*), which supports the left multiplication with a scalar.
        '''
        result=HuBbardList()
        for obj in self:
            result.append(obj*other)
        return result

    def __rmul__(self,other):
        '''
        Overloaded operator(*), which supports the right multiplication with a scalar.
        '''
        return self.__mul__(other)

    def bond_2nd_order(self,bonds):
        result=Superbonddict()
        for term in self:
            result[term.tag]=Superbondlist()
        for bond in bonds:            
            for term in self:
                if bond.n==2:
                    if hasattr(term,'point'):                     
                        if callable(term.point):
                            if term.point(bond.points[0]):
                                result[term.tag].append(bond)
                            if term.point(bond.points[1]):
                                result[term.tag].append(bond.reverse)
                        else:
                            if bond.sites[0] in term.point:
                                result[term.tag].append(bond)
                            if bond.sites[1] in term.point:
                                result[term.tag].append(bond.reverse)
                    else:
                        result[term.tag].append(bond)
                        result[term.tag].append(bond.reverse)                       
        return result

