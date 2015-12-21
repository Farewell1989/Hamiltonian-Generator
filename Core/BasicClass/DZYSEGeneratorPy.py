# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 19:35:20 2015

@author: DZY
"""
from numpy import *
from ConstantPy import RZERO
#from OperatorPy import *
from numpy.linalg import norm
from scipy.sparse import csr_matrix,coo_matrix
from scipy import sparse
from collections import OrderedDict

class SEGenerator:
    '''
    This class provides methods to generate and update operators of a Hamiltonian according to the bonds of the model and the descriptions of its terms.
    Attributes:
        bonds: list of Bond
            The bonds of the model.
        table: Table
            The index-sequence table of the system
        parameters: dict
            It contains all model parameters, divided into two groups, the constant ones and the alterable ones.
        terms: dict
            It contains all terms contained in the model, divided into two groups, the constant ones and the alterable ones.
        nambu: logical
            A flag to tag whether the nambu space takes part.
        half: logical
            A flag to tag whether the generated operators contains its hermitian conjugate part.
        cache: dict
            The working space used to handle the generation and update of the operators.
    '''

    def __init__(self,latticedim,bonds,terms=None):
        '''
        Constructor.
        Parameter:
            bonds: list of Bond
            latticedim: the Lattice.dim
            terms: list of Term,optional
                All the terms that are contained in the model.
                Those terms having the attribute modulate will go into self.terms['alter'] and the others will go into self.terms['const'].
            nambu: logical,optional
            half: logical,optional
        '''
        self.bonds=terms.bondclassify(bonds)
        self.latticedim=latticedim
        self.parameters={}
        self.terms={}
        self.set_parameters_and_terms(terms)
        self.cache={}
        self.set_cache()

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
        dim=prod(self.latticedim.values())
        if 'const' in self.terms:
            self.cache['const']=coo_matrix((dim,dim))
            for tag in self.terms['const'].iterkeys():
                for bond in self.bonds[tag]:
                    self.cache['const']=self.cache['const']+self.terms['const'][tag].hamiltonian(self.latticedim,bond)
        if 'alter' in self.terms:
            for tag in self.terms['alter'].iterkeys():
                self.cache['alter'][tag]=coo_matrix((dim,dim))
                for bond in self.bonds:
                    self.cache['alter'][tag]=self.cache['alter'][tag]+self.terms['alter'][tag].hamiltonian(self.latticedim,bond)

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

    @property
    def hamiltonian(self):
        '''
        This method returns the whole Hamiltonian of the systom.
        '''
        dim=prod(self.latticedim.values())
        result=coo_matrix((dim,dim))
        if 'const' in self.cache:
            result=result+self.cache['const']
        if 'alter' in self.cache:
            for opts in self.cache['alter'].itervalues():
                result=result+opts
        return result+result.conj().T

    def update(self,**karg):
        '''
        This method updates the alterable part of the Hamiltonian by keyword arguments.
        '''
        if 'alter' in self.terms:
            masks={key:1 for key in self.terms['alter'].iterkeys()}
            for key,term in self.terms['alter'].iteritems():
                nv=term[0].modulate(**karg)
                if nv is not None and norm(array(nv)-array(term[0].value))>RZERO:
                    masks[key]=nv/term[0].value
                    term[0].value=nv
                    self.parameters['alter'][key]=nv
            for key,mask in masks.iteritems():
                if mask!=1:
                    self.cache['alter'][key]=mask*self.cache['alter'][key]