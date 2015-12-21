# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 13:28:03 2015

@author: DZY
"""


from numpy import *
from copy import deepcopy

class SindexPack:
    def __init__(self,value,n=0,dim=[],index=[]):
        self.value=value
        self.n=n
        self.dim=dim
        self.index=index
        
    def __add__(self,other):
        result=SindexPackList()
        result.append(deepcopy(self))
        if isinstance(other,SindexPack):
            result.append(deepcopy(other))
        elif isinstance(other,SindexPackList):
            result.extend(deepcopy(other))
        else:
            raise ValueError("SindexPack error: the 'other' parameter must be of class SindexPack or SindexPackList.")
        return result
    
    def _mul(self,other):
        if isinstance(other,SindexPack):
            result=SindexPack(self.value*other.value,self.n+other.n,self.dim+other.dim,self.index+other.index)      
        else:
            result=deepcopy(self)
            result.value=self.value*other
        return SindexPackList(result)
    
    def __mul__(self,other):
        result=SindexPackList()
        if isinstance(other,SindexPack) or type(other)==int or type(other)==long or type(other)==float or type(other)==complex:
            result.extend(self._mul(other))
        elif isinstance(other,SindexPackList):
            for buff in other:
                result.extend(self._mul(buff))
        else:
            raise ValueError("SindexPack '*' error: the 'other' parameter must be of class SindexPack or .")
        return result
    
    def tomat(self):
        result=1.0
        for x in xrange(self.n):
            m=zeros((self.dim[x],self.dim[x]))
            m[self.index[x]]=1.0
            result=kron(result,m)
        return self.value*(result+result.conj().T)
        
class SindexPackList(list):
    def __init__(self,*arg):
        for buff in arg:
            if isinstance(buff,SindexPack):
                self.append(buff)
            else:
                raise ValueError("SindexPackList init error: the input parameters must be of class SindexPack.")
                
    def __add__(self,other):
        result=SindexPackList(*deepcopy(self))
        if isinstance(other,SindexPack):
            result.append(deepcopy(other))
        elif isinstance(other,SindexPackList):
            result.extend(deepcopy(other))
        else:
            raise ValueError("SindexPackList '+' error: the 'other' parameter must be of class SindexPack or SindexPackList.")
        return result
    
    def __mul__(self,other):
        result=SindexPackList()
        if isinstance(other,SindexPack) or isinstance(other,SindexPackList) or type(other)==int or type(other)==long or type(other)==float or type(other)==complex:
            for buff in self:
                result.extend(buff*other)
        else:
            raise ValueError("SindexPackList '*' error: the 'other' parameter must be of class SindexPack or SindexPackList.")
        return result
        
    def tomat(self):
        result=self[0].tomat()
        for buff in self[1:]:
            result=result+buff.tomat()
        return result

def Gmma(n,index):
    result=SindexPack(1.0,1,[n],[index])
    return SindexPackList(result)

