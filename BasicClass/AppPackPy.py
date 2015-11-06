from EngineAppPy import *
from numpy import *
class EB(App):
    '''
    Energy bands.
    '''
    def __init__(self,path=None,**karg):
        self.path=path
        self.emax=10.0 if 'emax' not in karg else karg['emax']
        self.emin=-10.0 if 'emin' not in karg else karg['emin']
        self.ne=400 if 'ne' not in karg else karg['ne']
        self.delta=0.05 if 'delta' not in karg else karg['delta']
        if 'ts' in karg: self.ts=karg['ts']

class DOS(App):
    '''
    Density of states.
    '''
    def __init__(self,BZ=None,ne=100,delta=0.05,emin=-10.0,emax=10.0,**karg):
        self.BZ=BZ
        self.ne=ne
        self.delta=delta
        self.emin=emin
        self.emax=emax

class GP(App):
    '''
    Grand potential.
    '''
    def __init__(self,BZ=None,e1=5.0,e2=50.0,deg1=64,deg2=64,deg3=64,n=20,**karg):
        self.BZ=BZ
        self.e_degs=[(0,e1,deg1),(e1,e2,deg2),(e2,e2*n,deg3)]
        self.gp=0

class GPS(App):
    '''
    Grand potential surface.
    '''
    def __init__(self):
        pass

class GFC(App):
    '''
    The coefficients of Green's functions.
    '''
    def __init__(self,nstep=200,vtype='rd',**karg):
        self.nstep=nstep
        self.vtype=vtype
        self.gse=0
        self.coeff=array([])

class GF(App):
    '''
    Green's functions.
    '''
    def __init__(self,shape,omega=0.0,k=[],**karg):
        self.omega=omega
        self.k=k
        self.gf=zeros(shape,dtype=complex128)
