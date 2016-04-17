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
from DZYSpinWmatrixPy import *
#from TablePy import *
from QuadraticPy import * 
from ConstantPy import RZERO
from copy import deepcopy

class Superexchange(Term):
    def __init__(self,mode,tag,value,bond,mat=None,sindexpack=None,amplitude=None,modulate=None):
        super(Superexchange,self).__init__(mode,tag,value,modulate)
        self.amplitude=amplitude
        self.bond=bond
        dim=bond['dim']
        self.mat=zeros((dim,dim))
        self.sindexpack=None
        if mat is not None:
            if all(abs(mat-mat.conj().T)<=RZERO):
                self.mat=self.mat+(mat+mat.conj().T)/2.0
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

class SOC(Superexchange):
    def __init__(self,mode,tag,value,bond,mat=None,sindexpack=None,amplitude=None,modulate=None):
        if bond['test']!=0:
            raise ValueError('SOC error: the SOC term must define on only one site')
        super(SOC,self).__init__(mode,tag,value,bond,mat,sindexpack,amplitude,modulate)
        
    def to_GUterm(self,table):
        table=reverse_table(table)
        mat=self.mat
        seq=transpose(where(abs(mat)>0.1))
        indexpackages=IndexPackageList()
        for x in seq:
            i,j=x
            indexpackages.append(IndexPackage(value=mat[i,j],orbital1=table[i].orbital,orbital2=table[j].orbital,spin1=table[i].spin,spin2=table[j].spin))
        return Quadratic(self.mode,self.tag,self.value,neighbour=0,indexpackages=indexpackages,amplitude=None,modulate=None)
        

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
        dim=bond.dim
        mat=self.mesh(bond).reshape(dim+dim)
        x=where(mat!=0)
        h=_hamiltonian(latticedim,bond.series,mat[x],x)
        return h
        
    def spinWmat(self,bond):
        return SpinWmatrix(bond,deepcopy(self.mesh(bond)))
        
    def bondclassify(self,bonds):
        result=Superbonddict()
        for term in self:
            result[term.tag]=Superbondlist()
        for bond in bonds:            
            for term in self:
                flag=False
                if callable(term.bond['test']):
                    flag=term.bond['test'](bond)
                else:
                    flag=term.bond['test']==bond.tag
                if flag:
                    result[term.tag].append(bond)
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
    

def HUtoSE(t,ci,cj,p0i,p0j,p1i,p1j,mode,tag,bond,reverse=False):
    '''
    t: narray 2-d matrix
    hopping term. 
    ci,cj: dict{seq:narray matrix}
    annihilation operator
    p0i,p0j: dict{E:[phi1,phi2,...]}
    low energy projector
    p0i,p0j: dict{E:[phi1,phi2,...]}
    exciting energy projector
    '''
    si={}
    for e1i,phi1i in p1i.items():
        termi={}
        for e0i,phi0i in p0i.items():
            termi.update({e0i:_cmat(phi1i,phi0i,ci)})
        si.update({e1i:termi})
    sj={}
    for e1j,phi1j in p1j.items():
        termj={}
        for e0j,phi0j in p0j.items():
            termj.update({e0j:_cmat(phi0j,phi1j,cj)})
        sj.update({e1j:termj})                       
    return HUtoSEsimple(t,si,sj,mode,tag,bond,reverse)

def _cmat(phi0,phi1,c):
    out={}
    for x,ci in c.iteritems():
        mat=None
        for p0 in phi0:
            for p1 in phi1:
                if mat is None:
                    mat=(p0.conj().T.dot(ci).dot(p1))*p0.dot(p1.conj().T)
                else:
                    mat=mat+(p0.conj().T.dot(ci).dot(p1))*p0.dot(p1.conj().T)
        out.update({x:mat})
    return out
    
#def _genmat(phi0i,phi1i,phi0j,phi1j,cj,cjd,t,cid,ci,td):
#    out=None
#    for p0i in phi0i:
#        for p1i in phi1i:
#            for p0j in phi0j:
#                for p1j in phi1j:
#                    if out is None:
#                        out=t*td*(p0i.conj().T.dot(cid).dot(p1i))*(p1i.conj().T.dot(ci).dot(p0i))*(p0j.conj().T.dot(cj).dot(p1j))*(p1j.conj().T.dot(cjd).dot(p0j))*kron(p0i.dot(p0i.conj().T),p0j.dot(p0j.conj().T))
#                    else:
#                        out=out+t*td*(p0i.conj().T.dot(cid).dot(p1i))*(p1i.conj().T.dot(ci).dot(p0i))*(p0j.conj().T.dot(cj).dot(p1j))*(p1j.conj().T.dot(cjd).dot(p0j))*kron(p0i.dot(p0i.conj().T),p0j.dot(p0j.conj().T))
#    return out

def ci_cj(ci,cj,p0i,p0j,p1i,p1j):
    '''
    ci,cj: dict{seq:narray matrix}
    annihilation operator
    p0i,p0j: dict{E:[phi1,phi2,...]}
    low energy projector
    p0i,p0j: dict{E:[phi1,phi2,...]}
    exciting energy projector
    '''
    si={}
    for e1i,phi1i in p1i.items():
        termi={}
        for e0i,phi0i in p0i.items():
            termi.update({e0i:_cmat(phi1i,phi0i,ci)})
        si.update({e1i:termi})
    sj={}
    for e1j,phi1j in p1j.items():
        termj={}
        for e0j,phi0j in p0j.items():
            termj.update({e0j:_cmat(phi0j,phi1j,cj)})
        sj.update({e1j:termj})                       
    return si,sj

def HUtoSEsimple(t,si,sj,mode,tag,bond,reverse=False):
    '''
    t: narray 2-d matrix
    hopping term. 
    termi,termj: dict((E0,E1):{seq:narray matrix})
    annihilation operator
    '''
    mat=None
    for e1i,termi in si.iteritems():
        for e1j,termj in sj.iteritems():
            for eil,cid in termi.items():
                for eir,ci in termi.items():
                    for ejl,cj in termj.items():
                        for ejr,cjd in termj.items():
                            value=0.5/(eil+ejl-e1i-e1j)+0.5/(eir+ejr-e1i-e1j)
                            i=len(cid)
                            j=len(cj)
                            ii=len(ci)
                            jj=len(cjd)
                            mat1=None
                            for x in xrange(i):
                                for xx in xrange(ii):
                                    for y in xrange(j):
                                        for yy in xrange(jj):
                                            if mat1 is None:
                                                mat1=t[x,y]*t.conj().T[yy,xx]*kron(ci[x].conj().T.dot(ci[xx]),cj[y].dot(cj[yy].conj().T))
                                            else:
                                                mat1=mat1+t[x,y]*t.conj().T[yy,xx]*kron(ci[x].conj().T.dot(ci[xx]),cj[y].dot(cj[yy].conj().T))
                                            if reverse:
                                                mat1=mat1+t[yy,xx]*t.conj().T[x,y]*kron(cj[y].dot(cj[yy].conj().T),ci[x].conj().T.dot(ci[xx]))
                            if mat is None:
                                mat=value*mat1
                            else:
                                mat=mat+value*mat1
    return Superexchange(mode,tag,value=1.0,bond=bond,mat=mat)
    