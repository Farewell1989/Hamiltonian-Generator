# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 13:24:20 2015

@author: DZY
"""
from BasicGeometryPy import *
from copy import deepcopy
from collections import OrderedDict

class Superbond:
    def __init__(self,tag,points):
        self.tag=tag
        self.points=points
        self.n=len(points)
        self.sites=[point.site for point in points]
        
    @property
    def number(self):
        return self.n

    @property
    def series(self):
        temp={point.site:x for x,point in enumerate(self.points)}
        return OrderedDict(sorted(temp.items(), key=lambda x:x[1]))
        
    @property
    def dim(self):
        return [point.struct.dim() for point in self.points]

    @property
    def rcoords(self):
        return [point.rcoord for point in self.points]

    @property
    def icoords(self):
        return [point.icoord for point in self.points]

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