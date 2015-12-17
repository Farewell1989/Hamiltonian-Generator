'''
App pack.
'''
from EngineAppPy import *
from numpy import *
class EB(App):
    '''
    Energy bands.
    '''
    def __init__(self,path=None,**karg):
        '''
        Constructor.
        Parameters:
            path: BaseSpace, optional
                The path in basespace along which the energy spectrum is to be computed.
        '''
        self.path=path
        self.emax=10.0 if 'emax' not in karg else karg['emax']
        self.emin=-10.0 if 'emin' not in karg else karg['emin']
        self.ne=400 if 'ne' not in karg else karg['ne']
        self.eta=0.05 if 'eta' not in karg else karg['eta']
        self.ns=6 if 'ns' not in karg else karg['ns']
        if 'ts' in karg: self.ts=karg['ts']

class DOS(App):
    '''
    Density of states.
    '''
    def __init__(self,BZ=None,ne=100,eta=0.05,emin=-10.0,emax=10.0,**karg):
        '''
        Constructor.
        Parameters:
            BZ: BaseSpace,optional
                The Brillouin zone upon which the energy levels are to be computed.
            emin,emax: float, optional
                They define the range of the energy within which the DOS is to be computed.
            ne: int, optional
                The number of sample points in the energy range defined by emin and emax.
            eta: float, optional
                The damping factor.
        '''
        self.BZ=BZ
        self.ne=ne
        self.eta=eta
        self.emin=emin
        self.emax=emax

class OP(App):
    '''
    Order parameter.
    '''
    def __init__(self,terms,BZ=None,p=1.0,**karg):
        '''
        Constructor.
        Parameters:
            term: Term
                The term representing the order parameter.
            BZ: BaseSpace, optional
                The Brillouin zone.
        '''
        self.terms=terms
        self.BZ=BZ
        self.ms=0
        self.ops=0
        self.p=p

    def matrix(self,bonds,table,nambu):
        '''
        '''
        pass

class FF(App):
    '''
    Filling factor.
    '''
    def __init__(self,BZ,p=1.0,**karg):
        self.BZ=BZ
        self.p=p
        self.filling=0

class CP(App):
    '''
    Chemical potential.
    '''
    def __init__(self,BZ=None,eta=0.05,error=10**-6,a=-20,b=20,p=1.0,cut=-10,deg1=200,deg2=200,deg3=200,n1=10**2,n2=10**5,**karg):
        self.BZ=BZ
        self.p=p
        self.eta=eta
        self.error=10**-6
        self.a=a
        self.b=b
        self.e_degs=[(n2*cut,n1*cut,deg1),(n1*cut,cut,deg2),(cut,0,deg3)]
        self.mu=0

class FS(App):
    '''
    Fermi surface.
    '''
    def __init__(self,BZ,eta=0.05,**karg):
        self.BZ=BZ
        self.eta=eta

class GP(App):
    '''
    Grand potential.
    '''
    def __init__(self,BZ=None,**karg):
        self.BZ=BZ
        self.gp=0

class GPS(App):
    '''
    Grand potential surface.
    '''
    def __init__(self,BS,**karg):
        self.BS=BS

class CN(App):
    '''
    Chern number.
    '''
    def __init__(self,BZ,d=10**-6,**karg):
        self.BZ=BZ
        self.d=d
        self.bc=None

    @property
    def cn(self):
        return sum(self.bc)*self.BZ.volume['k']/self.BZ.rank['k']/2/pi

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
