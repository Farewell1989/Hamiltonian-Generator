# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 13:24:20 2015

@author: DZY
"""
from BasicGeometryPy import *
from copy import deepcopy

class Superbond:
    def __init__(self,tag,points):
        self.tag=tag
        self.n=len(points)
        self.sites={points[x].site:x for x in xrange(self.n)}
        self.points=points
    
    @property
    def number(self):
        return self.n
    
    @property
    def dim(self):
        result={}
        for x in self.points:
            result[x.site]=x.struct.dim()
        return result

    @property
    def reverse(self):
        points=deepcopy(self.points)
        points.reverse()
        return Superbond(self.tag,points)
    
    def sort(self,methed):
        points=sorted(self.points,key=methed)
        return Superbond(self.tag,points)


class Superbondlist(list):
    pass 

class Superbonddict(dict):
    pass