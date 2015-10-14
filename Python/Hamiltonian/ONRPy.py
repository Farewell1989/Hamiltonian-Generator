from BasicClass.AppPackPy import *
from BasicClass.QuadraticPy import *
from BasicClass.HubbardPy import *
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
    1) ensemble: 'c' for canonical ensemble and 'g' for grand canonical ensemble;
    2) name: the name of system;
    3) lattice: the lattice of the system;
    4) filling: the filling factor of the system;
    5) mu: the chemical potential of the system;
    6) basis: the occupation number basis of the system;
    7) nspin: a flag to tag whether the ground state of the system lives in the subspace where the spin up electrons equal the spin down electrons, 1 for yes and 2 for no; 
    8) hopping: the hopping terms of the system, which is an instance of QuadraticList;
    9) onsite: the onsite terms of the system, which is an instance of QuadraticList;
    10) pairing: the pairing terms of the system, which is an instance of QuadraticList;
    11) hubbard: the Hubbard interaction terms of the system, which is an instance of HubbardList;
    12) boundary: the boundary condition of the system, 'op' for open and 'pr' for periodic;
    13) nambu: a flag to tag whether the system works in the nambu space;
    14) table: the index-sequence table of the system;
    15) operators: a dict containing different groups of operators for diverse tasks, e.g. entry 'h' includes "half" the operators of the Hamiltonian, entry 'gf' includes all the single-fermion operators each of which represents the index of the Green function, etc;
    16) matrix: the sparse matrix representation of the system.
    '''

    def __init__(self,ensemble='c',name=None,lattice=None,filling=0.5,mu=0,basis=None,nspin=1,hopping=[],onsite=[],pairing=[],hubbard=[],boundary='op',nambu=False,**karg):
        self.ensemble=ensemble
        self.name=Name(name)
        self.lattice=lattice
        self.filling=filling
        self.mu=mu
        self.basis=basis
        if basis.basis_type=='ES':
            self.nspin=nspin
        else:
            self.nspin=2
        self.hopping=HoppingList(*hopping)
        self.onsite=OnsiteList(*onsite)
        self.pairing=PairingList(*pairing)
        self.hubbard=HubbardList(*hubbard)
        self.boundary=boundary
        self.nambu=nambu
        self.table=Table(self.lattice.indices(nambu=self.nambu))
        self.operators={}

    def get_ready(self):
        self.set_name()
        self.set_operators()

    def set_name(self):
        self.name=Name(prefix=self.name.prefix,suffix=self.__class__.__name__)
        if self.ensemble.lower()=='c':
            self.name.append(self.filling)
        elif self.ensemble.lower()=='g':
            self.name.append(self.mu)
        self.name.extend(self.hopping)
        self.name.extend(self.onsite)
        self.name.extend(self.pairing)
        self.name.extend(self.hubbard)

    def set_operators(self):
        self._set_operators_hamiltonian()
        self._set_operators_greenfunction()

    def _set_operators_hamiltonian(self):
        nhp=len(self.hopping)
        nst=len(self.onsite)
        npr=len(self.pairing)
        nhb=len(self.hubbard)
        self.operators['h']=OperatorList()
        for neighbour,bonds in enumerate(self.lattice.bonds):
            for bond in bonds:
                if self.bond_fit_boundary(bond):
                    if neighbour==0:
                        if nhb>0: self.operators['h'].extend(self.hubbard.operators(bond,self.table,half=True))
                        if nst>0: self.operators['h'].extend(self.onsite.operators(bond,self.table,half=True))
                    else:
                        if nhp>0: self.operators['h'].extend(self.hopping.operators(bond,self.table,half=True))
                        if npr>0: self.operators['h'].extend(self.pairing.operators(bond,self.table,half=True))

    def _set_operators_greenfunction(self):
        self.operators['gf']=OperatorList()
        table=self.table if self.nspin==2 else subset(self.table,mask=lambda index: True if index.spin==0 else False)
        for index,sequence in table.iteritems():
            if isinstance(index,Index):self.operators['gf'].append(E_Linear(1,indices=[index],rcoords=[self.lattice.points[index.site].rcoord],icoords=[self.lattice.points[index.site].icoord],seqs=[sequence]))
        self.operators['gf'].sort(key=lambda operator: operator.seqs[0])

    def bond_fit_boundary(self,bond):
        if self.boundary.lower()=='op':
            return bond.is_intra_cell()
        elif self.boundary.lower()=='pr':
            return True
        else:
            raise ValueError("ONR bond_fit_boundary error: boundary type '"+self.boundary+"' is not supported.")

    def set_matrix(self):
        self.matrix=csr_matrix((self.basis.nbasis,self.basis.nbasis),dtype=GP.S_dtype)
        for operator in self.operators['h']:
            self.matrix+=opt_rep(operator,self.basis,transpose=False)
        self.matrix+=conjugate(transpose(self.matrix))
        self.matrix=transpose(self.matrix)

def ONRGFC(engine,app):
    nopt=len(engine.operators['gf'])
    if os.path.isfile(engine.din+'/'+engine.name.full_name+'_coeff.dat'):
        with open(engine.din+'/'+engine.name.full_name+'_coeff.dat','rb') as fin:
            app.gse=fromfile(fin,count=1)
            app.coeff=fromfile(fin,dtype=GP.S_dtype)
        if len(app.coeff)==nopt*nopt*2*3*app.nstep:
            app.coeff=app.coeff.reshape((nopt,nopt,2,3,app.nstep))
            return
    app.coeff=zeros((nopt,nopt,2,3,app.nstep),dtype=GP.S_dtype)
    engine.set_matrix()
    app.gse,gs=Lanczos(engine.matrix,vtype=app.vtype).eig(job='v')
    print 'gse:',app.gse
    if engine.basis.basis_type.lower() in ('es','ep'): engine.matrix=None
    for h in xrange(2):
        if h==0: print 'Electron part:'
        else: print 'Hole part:' 
        for j,optb in enumerate(engine.operators['gf']):
            for i,opta in enumerate(engine.operators['gf']):
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
        result.matrix=csr_matrix((result.basis.nbasis,result.basis.nbasis),dtype=GP.S_dtype)
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
        result.matrix=csr_matrix((result.basis.nbasis,result.basis.nbasis),dtype=GP.S_dtype)
        result.set_matrix()
        return result

def ONRGF(engine,app):
    nmatrix=engine.apps['GFC'].nstep
    gse=engine.apps['GFC'].gse
    coeff=engine.apps['GFC'].coeff
    nopt=len(engine.operators['gf'])
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
    engine.addapps(app=GF((len(engine.operators['gf']),len(engine.operators['gf'])),run=ONRGF))
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

# The following codes are used for tests only.
def test_onr():
    U=2.0
    t=-1.0
    m=2;n=4
    p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],atom=0,norbital=1,nspin=2,nnambu=1,scope='WG'+str(m)+str(n))
    a1=array([1.0,0.0])
    a2=array([0.0,1.0])
    a=ONR(
            din=        'Results/Coeff',
            dout=       'Results',
            ensemble=   'c',
            name=       'WG'+str(m)+str(n),
            lattice=    Lattice(name='WG'+str(m)+str(n),points=[p1],translations=[(a1,m),(a2,n)]),
            filling=    0.5,
            mu=         U/2,
            basis=      BasisE(up=(m*n,m*n/2),down=(m*n,m*n/2)),
#            basis=      BasisE((2*m*n,m*n)),
            nspin=      2,
            hopping=    [Hopping(t,neighbour=1)],
            hubbard=    [Hubbard([U])],
            boundary=   'op'
            )
    a.get_ready()
    print a.operators['gf']
    a.addapps('GFC',GFC(nstep=200,save_data=False,vtype='SY',run=ONRGFC))
    a.addapps('DOS',DOS(emin=-5,emax=5,ne=401,delta=0.05,run=ONRDOS,show=False))
    a.runapps()