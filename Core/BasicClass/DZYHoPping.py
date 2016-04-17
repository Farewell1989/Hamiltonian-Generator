# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 20:32:31 2016

@author: DZY
"""
from numpy import *
from TermPy import *
from DZYSuperbondPy import *
from DZYSpinWmatrixPy import *
from Hamiltonian.Core.BasicAlgorithm.KroneckerPy import *

class HoPping(Term):
    '''
    ''' 
    def __init__(self,mode,tag,value,bond,tspin=None,torbital=None,t=None,amplitude=None,modulate=None):
        '''
        Constructor.
        Parameters:
            mode: string
                The type of the term.
            tag: string
                The tag specifying the term used for dictionary lookup.
            value: float or complex
                The overall coefficient of the term.
            bond: bondtag or callable
                Indicate the bond the hoping term on.
            amplitude: function
                It must return a float or complex and take an instance of Bond as its only argument.
            modulate: function
                It must return a float or complex and its arguments are unlimited.
        '''
        super(HoPping,self).__init__(mode,tag,value,modulate)
        self.bond=bond
        self.tspin=tspin
        self.torbital=torbital
        if t is not None:
            self.t=t
        else:
            self.t=kron(tspin,torbital)
        self.amplitude=amplitude

#    def __str__(self):
#        '''
#        Convert an instance to string.
#        '''
#        result='Mode & tag: '+self.mode+' & '+self.tag+'\n'+'Value: '+str(self.value)+'\n'
#        if isinstance(self.indexpackages,IndexPackageList):
#            result+=str(self.indexpackages)
#        if hasattr(self,'modulate'):
#            result+='Modulate function: '+str(self.modulate)+'\n'
#        return result

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a Quadratic instance with a Quadratic/QuadraticList instance.
        '''
        result=HoPpingList()
        result.append(deepcopy(self))
        if isinstance(other,HoPping):
            result.append(deepcopy(other))
        elif isinstance(other,HoPpingList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('HoPping "+" error: the other parameter must be an instance of HoPping or HoPpingList.')
        return result

    def __pos__(self):
        '''
        Overloaded operator(+), i.e. +self.
        '''
        result=HoPpingList()
        result.append(deepcopy(self))
        return result


class HoPpingList(list):
    def __init__(self,*arg):
        for obj in arg:
            if isinstance(obj,HoPping):
                self.append(obj)
            else:
                raise ValueError("HoPpingList init error: the input argument should be HoPping instances.")

    def __add__(self,other):
        result=deepcopy(self)
        if isinstance(other,HoPping):
            result.append(deepcopy(other))
        elif isinstance(other,HoPpingList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('HoPpingList "+" error: the other parameter must be an instance of HoPping or HoPpingList.')
        return result

    def __mul__(self,other):
        result=HoPpingList()
        for obj in self:
            result.append(obj*other)
        return result

    def __rmul__(self,other):
        return self.__mul__(other)
        
    def spinW(self,spinWM):
        bond=spinWM.bond
        im=bond.points[0].struct.dim()
        jm=bond.points[1].struct.dim()
        dim=(im,jm)
        t=zeros(dim)
        for term in self:
            flag=term.bond(bond)
            if flag==1:
                t=t+term.value*term.t
            elif flag==-1:
                t=t+(term.value*term.t).conj().T
        if all(t==zeros(dim)):        
            t=eye(*dim)
        spinWM.t=kron(eye(im,im),t)
        spinWM.tp=kron(eye(im,im),t.conj().T)
        spinWM.mat=spinWM.tp.dot(spinWM.ormat).dot(spinWM.t)
        
#    def bondclassify(self,bonds):
#        result=Superbonddict()
#        for term in self:
#            result[term.tag]=Superbondlist()
#        for bond in bonds:            
#            for term in self:
#                flag=False
#                if callable(term.bond['test']):
#                    flag=term.bond['test'](bond)
#                else:
#                    flag=term.bond['test']==bond.tag
#                if flag:
#                    result[term.tag].append(bond)
#        return result