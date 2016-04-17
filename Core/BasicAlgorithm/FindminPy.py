# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 21:51:39 2015

@author: DZY
"""

from numpy import *
from scipy.optimize import fmin,fmin_bfgs

def findmin(f,xmin,xmax,turn,args=(),xtol=0.0001,ftol=0.0001,maxiter=None,maxfun=None,disp=0,retall=0,callback=None):
    l=len(xmin)
    opts = {'xtol': xtol,
            'ftol': ftol,
            'maxiter': maxiter,
            'maxfun': maxfun,
            'disp': False,
            'full_output': True,
            'retall': retall}
    fval=f(xmax)
    x=xmax
    for i in xrange(turn):
        x0=random.rand(l)*(xmax-xmin)+xmin
        res=fmin(f, x0, args, callback=callback,**opts)
        if res[1]<=fval:
            fval=res[1]
            x=res[0]
    if disp:
        opts['disp']=True
        res=fmin(f, x, args, callback=callback,**opts)
        print res
    return x,fval
       
    


def findmin_bfgs(f,xmin, xmax,turn,fprime=None, args=(), gtol=1e-05, norm=inf, epsilon=1.4901161193847656e-08, maxiter=None, disp=0, retall=0, callback=None):
    l=len(xmin)
    opts = {'gtol': gtol,
            'norm': norm,
            'epsilon': epsilon,
            'disp': False,
            'full_output':True,
            'maxiter': maxiter,
            'retall': retall}
    fval=f(xmax)
    x=xmax
    for i in xrange(turn):
        x0=random.rand(l)*(xmax-xmin)+xmin   
        res=fmin_bfgs(f, x0,fprime,args, callback=callback,**opts)
        if res[1]<=fval:
            fval=res[1]
            x=res[0]
    if disp:
        opts['disp']=True
        res=fmin_bfgs(f, x,fprime, args, callback=callback,**opts)
        print res
    return x,fval