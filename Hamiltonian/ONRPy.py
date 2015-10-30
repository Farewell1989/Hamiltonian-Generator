from BasicClass.AppPackPy import *
from BasicClass.QuadraticPy import *
from BasicClass.HubbardPy import *
from BasicClass.GeneratorPy import *
from BasicClass.NamePy import *
from BasicClass.BasisEPy import *
from BasicClass.OperatorRepresentationPy import *
from BasicAlgorithm.LanczosPy import *
from scipy.sparse.linalg import eigsh
from scipy.linalg import solve_banded,solveh_banded
from copy import deepcopy
import matplotlib.pyplot as plt
import os.path,sys
class ONR(Engine):
    '''
    The class ONR provides the methods to get the sparse matrix representation on the occupation number basis of an electron system. Apart from those inherited from its parent class Engine, it has the following attributes:
    1) name: the name of system;
    2) ensemble: 'c' for canonical ensemble and 'g' for grand canonical ensemble;
    3) filling: the filling factor of the system;
    4) mu: the chemical potential of the system;
    5) basis: the occupation number basis of the system;
    6) nspin: a flag to tag whether the ground state of the system lives in the subspace where the spin up electrons equal the spin down electrons, 1 for yes and 2 for no; 
    7) generator: the operator generator;
    8) operators: a dict containing different groups of operators for diverse tasks, e.g. entry 'h' includes "half" the operators of the Hamiltonian, entry 'sp' includes all the single-particle operators, etc;
    9) matrix: the sparse matrix representation of the system.
    '''

    def __init__(self,name=None,ensemble='c',filling=0.5,mu=0,basis=None,nspin=1,generator=None,**karg):
        self.name=Name(prefix=name,suffix=self.__class__.__name__)
        self.ensemble=ensemble
        self.filling=filling
        self.mu=mu
        if self.ensemble.lower()=='c':
            self.name['filling']=self.filling
        elif self.ensemble.lower()=='g':
            self.name['mu']=self.mu
        self.basis=basis
        if basis.basis_type=='ES':
            self.nspin=nspin
        else:
            self.nspin=2
        self.generator=generator
        self.name.update(self.generator.parameters['const'])
        self.name.update(self.generator.parameters['alter'])
        self.operators={}
        self.set_operators()

    def set_operators(self):
        self.set_operators_hamiltonian()
        self.set_operators_single_particle()

    def set_operators_hamiltonian(self):
        self.operators['h']=self.generator.operators

    def set_operators_single_particle(self):
        self.operators['sp']=OperatorList()
        table=self.generator.table if self.nspin==2 else subset(self.generator.table,mask=lambda index: True if index.spin==0 else False)
        for index,sequence in table.iteritems():
            if isinstance(index,Index):self.operators['sp'].append(E_Linear(1,indices=[index],rcoords=[self.generator.lattice.points[index.site].rcoord],icoords=[self.generator.lattice.points[index.site].icoord],seqs=[sequence]))
        self.operators['sp'].sort(key=lambda operator: operator.seqs[0])

    def set_matrix(self):
        self.matrix=csr_matrix((self.basis.nbasis,self.basis.nbasis),dtype=complex128)
        for operator in self.operators['h']:
            self.matrix+=opt_rep(operator,self.basis,transpose=False)
        self.matrix+=conjugate(transpose(self.matrix))
        self.matrix=transpose(self.matrix)

def ONRGFC(engine,app):
    nopt=len(engine.operators['sp'])
    if os.path.isfile(engine.din+'/'+engine.name.full_name+'_coeff.dat'):
        with open(engine.din+'/'+engine.name.full_name+'_coeff.dat','rb') as fin:
            app.gse=fromfile(fin,count=1)
            app.coeff=fromfile(fin,dtype=complex128)
        if len(app.coeff)==nopt*nopt*2*3*app.nstep:
            app.coeff=app.coeff.reshape((nopt,nopt,2,3,app.nstep))
            return
    app.coeff=zeros((nopt,nopt,2,3,app.nstep),dtype=complex128)
    engine.set_matrix()
    app.gse,gs=Lanczos(engine.matrix,vtype=app.vtype).eig(job='v')
    print 'gse:',app.gse
    if engine.basis.basis_type.lower() in ('es','ep'): engine.matrix=None
    for h in xrange(2):
        if h==0: print 'Electron part:'
        else: print 'Hole part:' 
        for j,optb in enumerate(engine.operators['sp']):
            for i,opta in enumerate(engine.operators['sp']):
                if engine.basis.basis_type.lower()=='es' and engine.nspin==2 and optb.indices[0].spin!=opta.indices[0].spin : continue
                mask=False
                if j==i and j==0 : mask=True
                if engine.basis.basis_type.lower()=='es' and engine.nspin==2 and j==i and j==nopt/2: mask=True
                if h==0:
                    if mask: onr=onr_eh(engine,optb.indices[0].dagger)
                    matj=opt_rep(optb.dagger,[engine.basis,onr.basis],transpose=True)
                    mati=opt_rep(opta.dagger,[engine.basis,onr.basis],transpose=True)
                else:
                    if mask: onr=onr_eh(engine,optb.indices[0])
                    matj=opt_rep(opta,[engine.basis,onr.basis],transpose=True)
                    mati=opt_rep(optb,[engine.basis,onr.basis],transpose=True)
                statei=mati.dot(gs)
                statej=matj.dot(gs)
                normj=norm(statej)
                statej[:]=statej[:]/normj
                lcz=Lanczos(onr.matrix,statej)
                for k in xrange(app.nstep):
                    if not lcz.cut:
                        app.coeff[i,j,h,0,k]=vdot(statei,statej)*normj
                        lcz.iter()
                app.coeff[i,j,h,1,0:len(lcz.a)]=array(lcz.a)
                app.coeff[i,j,h,2,0:len(lcz.b)]=array(lcz.b)
                print j*nopt+i,'...',
                sys.stdout.flush()
        print
    if app.save_data:
        with open(engine.din+'/'+engine.name.full_name+'_coeff.dat','wb') as fout:
            array(app.gse).tofile(fout)
            app.coeff.tofile(fout)

def onr_eh(self,index):
    if self.basis.basis_type.lower()=='eg':
        return self
    elif self.basis.basis_type.lower()=='ep':
        result=deepcopy(self)
        if index.nambu==CREATION:
            result.basis=BasisE((self.basis.nstate,self.basis.nparticle+1))
        else:
            result.basis=BasisE((self.basis.nstate,self.basis.nparticle-1))
        result.matrix=csr_matrix((result.basis.nbasis,result.basis.nbasis),dtype=complex128)
        result.set_matrix()
        return result
    else:
        result=deepcopy(self)
        if index.nambu==CREATION and index.spin==0:
            result.basis=BasisE(up=(self.basis.nstate[0],self.basis.nparticle[0]),down=(self.basis.nstate[1],self.basis.nparticle[1]+1))
        elif index.nambu==ANNIHILATION and index.spin==0:
            result.basis=BasisE(up=(self.basis.nstate[0],self.basis.nparticle[0]),down=(self.basis.nstate[1],self.basis.nparticle[1]-1))
        elif index.nambu==CREATION and index.spin==1:
            result.basis=BasisE(up=(self.basis.nstate[0],self.basis.nparticle[0]+1),down=(self.basis.nstate[1],self.basis.nparticle[1]))
        else:
            result.basis=BasisE(up=(self.basis.nstate[0],self.basis.nparticle[0]-1),down=(self.basis.nstate[1],self.basis.nparticle[1]))
        result.matrix=csr_matrix((result.basis.nbasis,result.basis.nbasis),dtype=complex128)
        result.set_matrix()
        return result

def ONRGF(engine,app):
    nmatrix=engine.apps['GFC'].nstep
    gse=engine.apps['GFC'].gse
    coeff=engine.apps['GFC'].coeff
    nopt=len(engine.operators['sp'])
    app.gf[...]=0.0
    buff=zeros((3,nmatrix),dtype=complex128)
    b=zeros(nmatrix,dtype=complex128)
    for i in xrange(nopt):
        for j in xrange(nopt):
            for h in xrange(2):
                b[...]=0;b[0]=1
                buff[0,1:]=coeff[i,j,h,2,0:nmatrix-1]*(-1)**(h+1)
                buff[1,:]=app.omega-(coeff[i,j,h,1,:]-gse)*(-1)**h
                buff[2,:]=coeff[i,j,h,2,:]*(-1)**(h+1)
                app.gf[i,j]+=inner(coeff[i,j,h,0,:],solve_banded((1,1),buff,b,overwrite_ab=True,overwrite_b=True,check_finite=False))

def ONRDOS(engine,app):
    engine.addapps(app=GF((len(engine.operators['sp']),len(engine.operators['sp'])),run=ONRGF))
    result=zeros((app.ne,2))
    for i,omega in enumerate(linspace(app.emin,app.emax,num=app.ne)):
        engine.apps['GF'].omega=omega+engine.mu+app.delta*1j
        engine.runapps('GF')
        result[i,0]=omega
        result[i,1]=-2*imag(trace(engine.apps['GF'].gf))
    if app.save_data:
        savetxt(engine.dout+'/'+engine.name.full_name+'_DOS.dat',result)
    if app.plot:
        plt.title(engine.name.full_name+'_DOS')
        plt.plot(result[:,0],result[:,1])
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full_name+'_DOS.png')
