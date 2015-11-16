from Hamiltonian.Core.BasicClass.AppPackPy import *
from Hamiltonian.Core.BasicClass.QuadraticPy import *
from Hamiltonian.Core.BasicClass.GeneratorPy import *
from Hamiltonian.Core.BasicClass.NamePy import *
from Hamiltonian.Core.BasicAlgorithm.BerryCurvaturePy import *
from scipy.linalg import eigh
import matplotlib.pyplot as plt 
class TBA(Engine):
    '''
    The TBA class provides a general algorithm to calculate physical quantities of non-interacting systems based on the tight-binding approximation. The BdG systems, i.e. phenomenological superconducting systems based on mean-field theory are also supported in a unified way. Apart from those inherited from its parent class Engine, TBA has the following attributes:
    1) filling: the filling factor of the system;
    2) mu: the chemical potential of the system;
    3) lattice: the lattice of the system;
    4) terms: the terms of the systems;
    5) generator: the operator generator;
    Supported apps and the corresponding run methods include:
    1) EB (energy bands) with method TBAEB,
    2) DOS (density of states) with TBADOS.
    '''
    
    def __init__(self,filling=0,mu=0,lattice=None,terms=None,nambu=False,**karg):
        self.filling=filling
        self.mu=mu
        self.lattice=lattice
        self.terms=terms
        self.generator=Generator(bonds=lattice.bonds,table=lattice.table(nambu),terms=terms,nambu=nambu,half=True)
        self.name.update(const=self.generator.parameters['const'])

    def matrix(self,k=[],**karg):
        self.generator.update(**karg)
        nmatrix=len(self.generator.table)
        result=zeros((nmatrix,nmatrix),dtype=complex128)
        for opt in self.generator.operators:
            phase=1 if len(k)==0 else exp(-1j*inner(k,opt.rcoords[0]))
            result[opt.seqs]+=opt.value*phase
            if self.generator.nambu:
                i,j=opt.seqs
                if i<nmatrix/2 and j<nmatrix/2: result[j+nmatrix/2,i+nmatrix/2]+=-opt.value*conjugate(phase)
        result+=conjugate(result.T)
        return result

    def eigvals(self,kspace=None):
        nmatrix=len(self.generator.table)
        result=zeros(nmatrix*(1 if kspace==None else product(kspace.rank.values())))
        if kspace==None:
            result[...]=eigh(self.matrix(),eigvals_only=True)
        else:
            for i,paras in enumerate(kspace()):
                result[i*nmatrix:(i+1)*nmatrix]=eigh(self.matrix(**paras),eigvals_only=True)
        return result

def TBAEB(engine,app):
    nmatrix=len(engine.generator.table)
    if app.path!=None:
        key=app.path.mesh.keys()[0]
        result=zeros((app.path.rank[key],nmatrix+1))
        if len(app.path.mesh[key].shape)==1:
            result[:,0]=app.path.mesh[key]
        else:
            result[:,0]=array(xrange(app.path.rank[key]))
        for i,parameter in enumerate(list(app.path.mesh[key])):
            result[i,1:]=eigh(engine.matrix(**{key:parameter}),eigvals_only=True)
    else:
        result=zeros((2,nmatrix+1))
        result[:,0]=array(xrange(2))
        result[0,1:]=eigh(engine.matrix(),eigvals_only=True)
        result[1,1:]=result[0,1:]
    if app.save_data:
        savetxt(engine.dout+'/'+engine.name.full+'_EB.dat',result)
    if app.plot:
        plt.title(engine.name.full+'_EB')
        plt.plot(result[:,0],result[:,1:])
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full+'_EB.png')
        plt.close()

def TBADOS(engine,app):
    result=zeros((app.ne,2))
    eigvals=engine.eigvals(app.BZ)
    for i,v in enumerate(linspace(eigvals.min(),eigvals.max(),num=app.ne)):
       result[i,0]=v
       result[i,1]=sum(app.eta/((v-eigvals)**2+app.eta**2))
    if app.save_data:
        savetxt(engine.dout+'/'+engine.name.full+'_DOS.dat',result)
    if app.plot:
        plt.title(engine.name.full+'_DOS')
        plt.plot(result[:,0],result[:,1])
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full+'_DOS.png')
        plt.close()

def TBACP(engine,app):
    nelectron=int(round(engine.filling*app.BZ.rank['k']*len(engine.generator.table)))
    eigvals=sort((engine.eigvals(app.BZ)))
    app.mu=(eigvals[nelectron]+eigvals[nelectron-2])/2
    engine.mu=app.mu
    print 'mu:',app.mu

def TBACN(engine,app):
    H=lambda kx,ky: engine.matrix(k=[kx,ky])
    app.bc=zeros(app.BZ.rank['k'])
    for i,paras in enumerate(app.BZ()):
        app.bc[i]=berry_curvature(H,paras['k'][0],paras['k'][1],engine.mu,d=app.d)
    print 'Chern number(mu):',app.cn,'(',engine.mu,')'
    if app.save_data or app.plot:
        buff=zeros((app.BZ.rank['k'],3))
        buff[:,0:2]=app.BZ.mesh['k']
        buff[:,2]=app.bc
    if app.save_data:
        savetxt(engine.dout+'/'+engine.name.full+'_BC.dat',buff)
    if app.plot:
        nk=int(round(sqrt(app.BZ.rank['k'])))
        plt.title(engine.name.full+'_BC')
        plt.axis('equal')
        plt.colorbar(plt.pcolormesh(buff[:,0].reshape((nk,nk)),buff[:,1].reshape((nk,nk)),buff[:,2].reshape((nk,nk))))
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full+'_BC.png')
        plt.close()
