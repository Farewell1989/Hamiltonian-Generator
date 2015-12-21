# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:15:40 2015

Superexchange Term

@author: DZY
"""

from TermPy import *
from DZYSindexPackPy import *
#from IndexPy import *
#from BondPy import *
from DZYSuperbondPy import *
#from TablePy import *
#from OperatorPy import * 

class Superexchange(Term):
    def __init__(self,mode,tag,value,bond,mat=None,sindexpack=None,amplitude=None,modulate=None):
        super(Superexchange,self).__init__(mode,tag,value,modulate)
        self.amplitude=amplitude
        self.bond=bond
        dim=bond['dim']
        self.mat=zeros((dim,dim))
        self.sindexpack=None
        if mat is not None:
            if all(mat==mat.conj().T):
                self.mat=self.mat+mat
            else:
                raise ValueError("Superexchange error: input mat should be Hermitian")
        elif sindexpack is not None:
            if isinstance(sindexpack,SindexPackList):
                self.mat=self.mat+sindexpack.tomat()             
            elif callable(sindexpack):
                self.sindexpack=sindexpack
        else:
            raise ValueError("Superexchange error: inputmat and sindexpack should not be both None.")

    def __add__(self,other):
        result=SuperexchangeList(deepcopy(self))
        if isinstance(other,Superexchange):
            result.append(deepcopy(other))
        elif isinstance(other,SuperexchangeList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('Superexchange "+" error: the other parameter must be an instance of Superexchange or SuperexchangeList')
        return result

    def __pos__(self):
        result=SuperexchangeList(deepcopy(self))
        return result

    def mesh(self,bond):
        value=self.value*(1 if self.amplitude==None else self.amplitude(bond))
        if callable(self.sindexpack):
            self.mat=self.mat+self.sindexpack(bond)
        result=triu(self.mat,1)
        result=result+diag(diag(self.mat)/2.0)
        return value*result
        

class SuperexchangeList(list):
    def __init__(self,*arg):
        for obj in arg:
            if isinstance(obj,Superexchange):
                self.append(obj)
            else:
                raise ValueError("SuperexchangeList init error: the input argument should be Superexchange instances.")

    def __add__(self,other):
        result=deepcopy(self)
        if isinstance(other,Superexchange):
            result.append(deepcopy(other))
        elif isinstance(other,SuperexchangeList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('SuperexchangeList "+" error: the other parameter must be an instance of Superexchange or SuperexchangeList.')
        return result

    def __mul__(self,other):
        result=SuperexchangeList()
        for obj in self:
            result.append(obj*other)
        return result

    def __rmul__(self,other):
        return self.__mul__(other)

    def mesh(self,bond):
        result=self[0].mesh(bond)
        for obj in self[1:]:
            result+=obj.mesh(bond)
        return result
        
    def hamiltonian(self,latticedim,bond):
        dim=bond.dim.values()
        mat=self.mesh(bond).reshape(dim+dim)
        x=where(mat!=0)
        h=_hamiltonian(latticedim,bond.sites,mat[x],x)
        return h
    
    def bondclassify(self,bonds):
        result=Superbonddict()
        for term in self:
            result[term.tag]=Superbondlist()
        for bond in bonds:
            flag=False
            for term in self:
                if callable(term.bond['test']):
                    flag=term.bond.test(bond)
                else:
                    flag=term.bond['test']==bond.tag
                if flag:
                    result[term.tag].append(bond)
                    break
        return result
            
                    
            
            
            

#class SuperexchangeList(dict):
#    def __int__(self,*arg):
#        for x in arg:
#            if x.bond['tag'] in self.keys():
#                self[x.bond['tag']].append(x)
#            else:
#                self[x.bond['tag']]=+x
            


      
from scipy.sparse import csr_matrix,coo_matrix
from scipy import sparse

def _hamiltonian(dim,points,values,dex):
    n=len(values)
    if n==0:
        v=product(dim.values())
        return coo_matrix((v,v))
    m=len(dex)/2
    mat=[1.0 for x in xrange(n)]
    for (i,v) in dim.items():
        if i in points.keys():
            for x in xrange(n):
                ma=coo_matrix(([1],([dex[points[i]][x]],[dex[m+points[i]][x]])),shape=(v,v))
                mat[x]=sparse.kron(ma,mat[x])
        else:
            mat=[sparse.kron(sparse.eye(v),mat[x]) for x in xrange(n)]
    mat1=values[0]*mat[0]
    for x in xrange(1,n):
        mat1=mat1+values[x]*mat[x]       
    return mat1