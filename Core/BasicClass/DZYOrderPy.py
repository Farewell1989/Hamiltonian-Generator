# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 21:42:55 2015

@author: DZY
"""

from numpy import *
from scipy.linalg import expm
from collections import OrderedDict

class order:
    def __init__(self,latticedim):
        self.bottom=OrderedDict()
        self.top=OrderedDict()
        for key,val in latticedim.iteritems():
            self.bottom[key]={'dim':val}
    
    def set_order(self,site,target,mat=None):
        self.top.update({site:{'target':target,'mat':mat}})
        del self.bottom[site]
        
    def settle(self):
        self.min=[]
        self.max=[]
        for key in self.bottom.iterkeys():
            dim=self.bottom[key]['dim']-1
            self.min=self.min+[0 for i in xrange(2*dim)]
            self.max=self.max+[pi/2 for i in xrange(dim)]+[2*pi for i in xrange(dim)]
            self.bottom[key].update({'mat':testCPS(dim)})
        self.min=array(self.min)
        self.max=array(self.max)
        keys=self.top.keys()
        for key in keys:
            if self.top[key]['target'] in keys:
                self._tobottom(self.top[key])
            
    def _tobottom(self,top):
        top['mat']=top['mat'].dot(self.top[top['target']]['mat'])
        top['target']=self.top[top['target']]['target']
        if top['target'] in self.top.keys():
            self._tobottom(top)
                
    def rotate(self,x):
        out={}
        for key in self.bottom.iterkeys():
            dim=2*(self.bottom[key]['dim']-1)
            x0=x[:dim]
            x=x[dim:]
            out.update({key:self.bottom[key]['mat'](x0)})
        for key in self.top.iterkeys():
            out.update({key:self.top[key]['mat'].dot(out[self.top[key]['target']])})
        return out
    
    def orderparameter(self,x,order):
        out={}
        rotate=self.rotate(x)
        for key,unitary in rotate.iteritems():
            out.update({key:unitary[:,0].conj().dot(order).dot(unitary[:,0])})
        return out


def testCPS(n):
    A=[]
    for i in xrange(n):
        m=zeros((n+1,n+1),dtype=complex)
        m[0,i+1]=1
        m[i+1,0]=1
        A.append(m)
        m=zeros((n+1,n+1),dtype=complex)
        m[0,i+1]=-1j
        m[i+1,0]=1j
        A.append(m)
    def out(x):
        sita=x[:n]
        psi=x[n:]
        X=[sita[0] for i in xrange(2*n)]
        for i in xrange(1,n):
            X[2*i-2]=X[2*i-2]*cos(sita[i])
            X[2*i-1]=X[2*i-1]*cos(sita[i])
        for i in xrange(n):
            X[2*i]=X[2*i]*cos(psi[i])
            X[2*i+1]=X[2*i+1]*sin(psi[i])
            for j in xrange(i+1,n):
                X[2*j]=X[2*j]*sin(sita[i+1])
                X[2*j+1]=X[2*j+1]*sin(sita[i+1])
        M=zeros((n+1,n+1))
        for i in xrange(2*n):
            M=M+X[i]*A[i]
        return expm(1j*M)
    return out      
        