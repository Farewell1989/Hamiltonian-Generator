# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 13:39:14 2015

@author: DZY
"""

from numpy import *

class SpinWmatrix:
    def __init__(self,bond,mat):
        '''
        sites: list
        mat:narray in numpy
        '''
        self.bond=bond
        self.mat=mat
        
    def update(self,rotate):
        unitary=array([1.0])
        for x in self.bond.sites:
            unitary=kron(unitary,rotate[x])
        self.mat=unitary.conj().T.dot(self.mat).dot(unitary)
    
    def testupdate(self,rotate):
        unitary=array([1.0])
        for x in self.bond.sites:
            unitary=kron(unitary,rotate[x])
        return real(unitary[:,0].conj().dot(self.mat+self.mat.conj().T).dot(unitary[:,0]))       
        
class SpinWmatrixList(list):
    
    pass
