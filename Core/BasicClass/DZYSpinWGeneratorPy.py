# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 20:32:59 2015

@author: DZY
"""

from numpy import *
from ConstantPy import RZERO
from DZYOrderPy import order
from DZYSpinWmatrixPy import *
from OperatorPy import *
from TablePy import *
from numpy.linalg import norm
from scipy.sparse import csr_matrix,coo_matrix
from scipy import sparse
from collections import OrderedDict
from Hamiltonian.Core.BasicAlgorithm.FindminPy import *

class SpinWGenerator:
    '''
    This class provides methods to generate and update operators of a Hamiltonian according to the bonds of the model and the descriptions of its terms.
    Attributes:
        bonds: list of Bond
            The bonds of the model.
        parameters: dict
            It contains all model parameters, divided into two groups, the constant ones and the alterable ones.
        terms: dict
            It contains all terms contained in the model, divided into two groups, the constant ones and the alterable ones.
        cache: dict
            The working space used to handle the generation and update of the operators.
    '''

    def __init__(self,latticedim,mag,bonds,terms=None):
        '''
        Constructor.
        Parameter:
            bonds: list of Bond
            mag: Class order
            terms: list of Term,optional
                All the terms that are contained in the model.
                Those terms having the attribute modulate will go into self.terms['alter'] and the others will go into self.terms['const'].
        '''
        self.latticedim=latticedim
        self.table=_lineartable(self.latticedim)
        self.bonds=bonds
        self.mag=mag
        self.parameters={}
        self.terms={}
        self.set_parameters_and_terms(terms)
        self.cache={}
        self.set_cache()
        self.findCLSE()
        self.linearoperators()

    def set_parameters_and_terms(self,terms):
        self.parameters['const']=OrderedDict()
        self.parameters['alter']=OrderedDict()
        self.terms['const']={}
        self.terms['alter']={}
        if not terms is None:
            for term in terms:
                if hasattr(term,'modulate'):
                    self.parameters['alter'][term.tag]=term.value
                    self.terms['alter'][term.tag]=+term
                else:
                    self.parameters['const'][term.tag]=term.value
                    self.terms['const'][term.tag]=+term
 

    def set_cache(self):
        if 'const' in self.terms:
            self.cache['const']=SpinWmatrixList()
            for tag in self.terms['const'].iterkeys():
                for bond in self.bonds[tag]:
                    self.cache['const'].append(self.terms['const'][tag].spinWmat(bond))
        if 'alter' in self.terms:
            self.cache['alter']={}
            for tag in self.terms['alter'].iterkeys():
                self.cache['alter'][tag]=SpinWmatrixList()
                for bond in self.bonds[tag]:
                    self.cache['alter'][tag].append(self.terms['alter'][tag].spinWmat(bond))

    def findCLSE(self,xmin=None,xmax=None,turn=10,**opt):
        f=self._classicalE
        xl=xmin if xmin is not None else self.mag.min
        xu=xmax if xmax is not None else self.mag.max
        x,energy=findmin_bfgs(f,xl,xu,turn,**opt)
        self.x=x
        self.CLSenergy=energy
        self._updatecache(x)
        return x,energy
        
    def _classicalE(self,x):
        rotate=self.mag.rotate(x)
        energy=0.0
        if 'const' in self.terms:
            for SpinWM in self.cache['const']:
                energy=energy+SpinWM.testupdate(rotate)
        if 'alter' in self.terms:
            for cache in self.terms['alter'].itervalues():
                for SpinWM in cache:
                    energy=energy+SpinWM.testupdate(rotate)
        return energy
        
    def orderparameter(self,orderoperator):
        return self.mag.orderparameter(self.x,orderoperator)
    
    def _updatecache(self,x):
        rotate=self.mag.rotate(x)
        if 'const' in self.terms:
            for SpinWM in self.cache['const']:
                SpinWM.update(rotate)
        if 'alter' in self.terms:
            for cache in self.cache['alter'].itervalues():
                for SpinWM in cache:
                    SpinWM.update(rotate)        


#    def __str__(self):
#        '''
#        Convert an instance to string.
#        '''
#        result='Parameters:\n'
#        if len(self.parameters['const']>0):
#            result+='Constant part: '+str(self.parameters['const'])+'\n'
#        if len(self.parameters['alter']>0):
#            result+='Alterable part: '+str(self.parameters['alter'])+'\n'
#        for key,obj in self.terms['const'].iteritems():
#            result+=key+':\n'+str(obj)
#        for key,obj in self.terms['alter'].iteritems():
#            result+=key+':\n'+str(obj)
#        return result
    
    def linearoperators(self):
        '''
        This method returns the linear operators.
        '''
        out=OperatorList()
        if 'const' in self.cache:
            for spinWM in self.cache['const']:
                out=out+_linearoperators(spinWM,self.table)
        if 'alter' in self.cache:
            for opts in self.cache['alter'].itervalues():
                for spinWM in opts:
                    out=out+_linearoperators(spinWM,self.table)
        self.operators=out

        
    def update(self,**karg):
        '''
        This method updates the alterable part of the Hamiltonian by keyword arguments.
        '''
        if 'alter' in self.terms:
            masks={key:False for key in self.terms['alter'].iterkeys()}
            for key,term in self.terms['alter'].iteritems():
                nv=term[0].modulate(**karg)
                if nv is not None and norm(array(nv)-array(term[0].value))>RZERO:
                    masks[key]=True
                    term[0].value=nv
                    self.parameters['alter'][key]=nv
            for tag,mask in masks.iteritems():
                if mask:
                    self.cache['alter'][tag]=SpinWmatrixList()
                    for bond in self.bonds[tag]:
                        self.cache['alter'][tag].append(self.terms['alter'][tag].spinWmat(bond))
            if any(masks.values()):
                self.findCLSE()
                self.linearoperators()

def _linearoperators(spinWM,table):
    out=OperatorList()
    m=len(table)
    num=spinWM.bond.n
    sites=spinWM.bond.sites
    dim=spinWM.bond.dim
    rcoords=spinWM.bond.rcoords
    icoords=spinWM.bond.icoords
    mat=spinWM.mat.reshape(dim+dim)
    for x in xrange(num):
        inded=[0 for i in xrange(num)]
        indee=[0 for i in xrange(num)]
        if abs(mat[tuple(inded+indee)])>RZERO:
            for d in xrange(1,dim[x]):
                out.append(Operator('onsite',-mat[tuple(inded+indee)],'tbc',[rcoords[x]-rcoords[x]],[rcoords[x]-rcoords[x]],(table[(sites[x],d)],table[(sites[x],d)])))
        for d in xrange(1,dim[x]):
            for e in xrange(1,dim[x]):
                inded[x]=d
                indee[x]=e
                if abs(mat[tuple(inded+indee)])>RZERO:
                    out.append(Operator('onsite',mat[tuple(inded+indee)],'tbc',[rcoords[x]-rcoords[x]],[rcoords[x]-rcoords[x]],(table[(sites[x],d)],table[(sites[x],e)])))
    for x in xrange(num-1):
        for y in xrange(x+1,num):
            inded=[0 for i in xrange(num)]
            indee=[0 for i in xrange(num)]
            for d in xrange(1,dim[x]):
                for e in xrange(1,dim[y]): 
                    inded[x]=d
                    inded[y]=e
                    if abs(mat[tuple(inded+indee)])>RZERO:
                        out.append(Operator('linear',mat[tuple(inded+indee)],'tbc',[rcoords[x]-rcoords[y]],[icoords[x]-icoords[y]],(table[(sites[x],d)],table[(sites[y],e)]+m)))
            inded=[0 for i in xrange(num)]
            indee=[0 for i in xrange(num)]
            for d in xrange(1,dim[x]):
                for e in xrange(1,dim[y]): 
                    inded[x]=d
                    indee[y]=e
                    if abs(mat[tuple(inded+indee)])>RZERO:
                        out.append(Operator('linear',mat[tuple(inded+indee)],'tbc',[rcoords[x]-rcoords[y]],[icoords[x]-icoords[y]],(table[(sites[x],d)],table[(sites[y],e)])))
            inded=[0 for i in xrange(num)]
            indee=[0 for i in xrange(num)]
            for d in xrange(1,dim[x]):
                for e in xrange(1,dim[y]): 
                    inded[y]=e
                    indee[x]=d
                    if abs(mat[tuple(inded+indee)])>RZERO:
                        out.append(Operator('linear',mat[tuple(inded+indee)],'tbc',[rcoords[y]-rcoords[x]],[icoords[y]-icoords[x]],(table[(sites[y],e)],table[(sites[x],d)])))
            inded=[0 for i in xrange(num)]
            indee=[0 for i in xrange(num)]
            for d in xrange(1,dim[x]):
                for e in xrange(1,dim[y]): 
                    indee[x]=d
                    indee[y]=e
                    if abs(mat[tuple(inded+indee)])>RZERO:
                        out.append(Operator('linear',mat[tuple(inded+indee)],'tbc',[rcoords[x]-rcoords[y]],[icoords[x]-icoords[y]],(table[(sites[x],d)]+m,table[(sites[y],e)])))
    return out

#def _linearoperators(spinWM,m):
#    out=OperatorList()
#    num=spinWM.bond.n
#    sites=spinWM.bond.sites
#    dim=spinWM.bond.dim
#    rcoords=spinWM.bond.rcoords
#    icoords=spinWM.bond.icoords
#    mat=spinWM.mat.reshape(dim+dim)
#    for x in xrange(num):
#        index=[0 for i in xrange(num)]
#        out.append(Operator('onsite',-mat[tuple(index+index)]*eye(dim[x]-1),'tbc',[0.0],[0.0],(sites[x],sites[x])))
#        index[x]=slice(1,dim[x])
#        out=out+Operator('onsite',mat[tuple(index+index)],'tbc',[0.0],[0.0],(sites[x],sites[x]))
#    for x in xrange(num-1):
#        for y in xrange(x+1,num):
#            inde1=[0 for i in xrange(num)]
#            inde2=[0 for i in xrange(num)]
#            inde1[x]=slice(1,dim[x])
#            inde1[y]=slice(1,dim[y])
#            out.append(Operator('linear',mat[tuple(inde1+inde2)],'tbc',[rcoords[x]-rcoords[y]],[icoords[x]-icoords[y]],(sites[x],sites[y]+m)))
#            inde1=[0 for i in xrange(num)]
#            inde2=[0 for i in xrange(num)]
#            inde1[x]=slice(1,dim[x])
#            inde2[y]=slice(1,dim[y])
#            out.append(Operator('linear',mat[tuple(inde1+inde2)],'tbc',[rcoords[x]-rcoords[y]],[icoords[x]-icoords[y]],(sites[x],sites[y])))
#            inde1=[0 for i in xrange(num)]
#            inde2=[0 for i in xrange(num)]
#            inde2[x]=slice(1,dim[x])
#            inde1[y]=slice(1,dim[y])
#            out.append(Operator('linear',mat[tuple(inde1+inde2)],'tbc',[rcoords[y]-rcoords[x]],[icoords[y]-icoords[x]],(sites[y],sites[x])))
#            inde1=[0 for i in xrange(num)]
#            inde2=[0 for i in xrange(num)]
#            inde2[x]=slice(1,dim[x])
#            inde2[y]=slice(1,dim[y])
#            out.append(Operator('linear',mat[tuple(inde1+inde2)],'tbc',[rcoords[x]-rcoords[y]],[icoords[x]-icoords[y]],(sites[x]+m,sites[y])))
#    return out

def _lineartable(latticedim):
    out=[]
    for i,dim in latticedim.items():
        index=[]
        for x in xrange(1,dim):
            index.append((i,x))
        out.append(Table(indices=index))   
    return union(out)