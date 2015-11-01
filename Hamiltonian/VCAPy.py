from ONRPy import *
from BasicClass.BaseSpacePy import *
from VCA_Fortran import *
class VCA(ONR):
    '''
    The class VCA implements the algorithm of the variational cluster approximation of an electron system. Apart from those inherited from the class Engine, it has the following attributes:
    1) name: the name of the system;
    2) ensemble: 'c' for canonical ensemble and 'g' for grand canonical ensemble;
    3) filling: the filling factor of the system;
    4) mu: the chemical potential of the system;
    5) basis: the occupation number basis of the system;
    6) nspin: a flag to tag whether the ground state of the system lives in the subspace where the spin up electrons equal the spin down electrons, 1 for yes and 2 for no;
    7) cell : the unit cell of the system;
    8) lattice: the cluster of the system;
    9) terms: the terms of the system;
    10) weiss: the Weiss terms added to the system;
    11) nambu: a flag to tag whether pairing terms are involved;
    12) generators: a dict containing the needed operator generators, which generally has two entries:
        (1) entry 'h' is the generator for the cluster Hamiltonian including weiss terms;
        (2) entry 'pt' is the generator for the perturbation terms including weiss ones;
    13) operators: a dict containing different groups of operators for diverse tasks, which generally has four entries:
        (1) entry 'h' includes "half" the operators of the Hamiltonian intra the cluster,
        (2) entry 'pt' includes "half" the operators of the perturbation terms inter the clusters,
        (3) entry 'sp' includes all the single-particle operators intra the cluster, and
        (4) entry 'csp' includes all the single-particle operators intra the unit cell;
    14) clmap: a dict containing the information needed to restore the translation symmetry broken by the choosing of the clusters, which has two entries:
        (1) 'seqs': a two dimensinal array whose element[i,j] represents the index sequence of the j-th single-particle operator within the cluster which should correspond to the i-th single-particle operator within the unit cell after the restoration of the translation symmetry;
        (2) 'coords': a three dimensinal array whose element[i,j,:] represents the rcoords of the j-th single-particle operator within the cluster which should correspond to the i-th single-particle operator within the unit cell after the restoration of the translation symmetry;
    15) matrix: the sparse matrix representation of the system;
    16) cache: the cache during the process of calculation.
    '''
    def __init__(self,name=None,ensemble='c',filling=0.5,mu=0,basis=None,nspin=1,cell=None,lattice=None,terms=None,weiss=None,nambu=False,**karg):
        self.name=Name(prefix=name,suffix=self.__class__.__name__)
        self.ensemble=ensemble
        self.filling=filling
        self.mu=mu
        if self.ensemble.lower()=='c':
            self.name['filling']=self.filling
        elif self.ensemble.lower()=='g':
            self.name['mu']=self.mu
        self.basis=basis
        self.nspin=nspin if basis.basis_type=='ES' else 2
        self.cell=cell
        self.lattice=lattice
        self.terms=terms
        self.weiss=weiss
        self.nambu=nambu
        self.generators={}
        self.generators['h']=Generator(
                    bonds=      [bond for bond in lattice.bonds if bond.is_intra_cell()],
                    table=      Table(lattice.indices(nambu=False)),
                    terms=      terms if weiss is None else terms+weiss,
                    nambu=      False,
                    half=       True
                    )
        self.generators['pt']=Generator(
                    bonds=      [bond for bond in lattice.bonds if not bond.is_intra_cell()],
                    table=      Table(lattice.indices(nambu=nambu)),
                    terms=      terms if weiss is None else terms+[term*(-1) for term in weiss],
                    nambu=      nambu,
                    half=       True
                    )
        self.name.update(self.generators['h'].parameters['const'])
        self.name.update(self.generators['h'].parameters['alter'])
        self.operators={}
        self.set_operators()
        self.clmap={}
        self.set_clmap()
        self.cache={}

    def set_operators(self):
        '''
        Prepare the operators needed in future calculations.
        Generally, there are four entries in the dict "self.operators":
        1) 'h': stands for 'Hamiltonian', which contains half of the intra-cluster operators of the Hamiltonian;
        2) 'pt': stands for 'perturbation', which contains half of the inter-cluster operators as the perturbation terms;
        3) 'sp': stands for 'single particle', which contains all the allowed or needed single particle operators within the cluster. When self.nspin==1 and self.basis.basis_type=='es' (spin-conserved systems), only spin-down single particle operators are included;
        4) 'csp': stands for 'cell single particle', which contains all the allowed or needed single particle operators within the unit cell. When self.nspin==1 and self.basis.basis_type=='es' (spin-conserved systems), only spin-down single particle operators are included.
        '''
        self.set_operators_hamiltonian()
        self.set_operators_perturbation()
        self.set_operators_single_particle()
        self.set_operators_cell_single_particle()

    def set_operators_perturbation(self):
        self.operators['pt']=OperatorList()
        table=self.generator['pt'].table if self.nspin==2 else subset(self.generators['pt'].table,mask=lambda index: True if index.spin==0 else False)
        for opt in self.generators['pt'].operators:
            if opt.indices[1] in table: self.operators['pt'].append(opt)

    def set_operators_cell_single_particle(self):
        self.operators['csp']=OperatorList()
        table=Table(self.cell.indices(nambu=self.nambu)) if self.nspin==2 else subset(Table(self.cell.indices(nambu=self.nambu)),mask=lambda index: True if index.spin==0 else False)
        for index,seq in table.iteritems():
            if isinstance(index,Index): 
                self.operators['csp'].append(E_Linear(1,indices=[index],rcoords=[self.cell.points[index.site].rcoord],icoords=[self.cell.points[index.site].icoord],seqs=[seq]))
        self.operators['csp'].sort(key=lambda opt: opt.seqs[0])

    def set_clmap(self):
        '''
        self.clmap is a dict which contains the necessary information to restore the translation symmetry broken by the different treating of inter and intra cluster quadratics.
        '''
        nsp,ncsp,ndim=len(self.operators['sp']),len(self.operators['csp']),len(self.operators['csp'][0].rcoords[0])
        buff=[]
        for i in xrange(ncsp):
            buff.append(OperatorList())
        for optl in self.operators['sp']:
            for i,optc in enumerate(self.operators['csp']):
                if optl.indices[0].orbital==optc.indices[0].orbital and optl.indices[0].spin==optc.indices[0].spin and optl.indices[0].nambu==optc.indices[0].nambu and has_integer_solution(optl.rcoords[0]-optc.rcoords[0],self.cell.vectors):
                    buff[i].append(optl)
                    break
        self.clmap['seqs'],self.clmap['coords']=zeros((ncsp,nsp/ncsp),dtype=int64),zeros((ncsp,nsp/ncsp,ndim),dtype=float64)
        for i in xrange(ncsp):
            for j,optj in enumerate(buff[i]):
                self.clmap['seqs'][i,j],self.clmap['coords'][i,j,:]=optj.seqs[0]+1,optj.rcoords[0]

    def pt(self,ks):
        '''
        Returns the matrix form of the perturbation terms.
        '''
        ngf=len(self.operators['sp'])
        result=zeros((ngf,ngf),dtype=complex128)
        for opt in self.operators['pt']:
            result[opt.seqs]+=opt.value*(1 if len(ks)==0 else exp(-1j*inner(ks,opt.icoords[0])))
        return result+conjugate(result.T)

    def pt_mesh(self,kmesh):
        '''
        Returns the mesh of the perturbation terms.
        '''
        if 'pt_mesh' in self.cache:
            return self.cache['pt_mesh']
        else:
            result=zeros((kmesh.shape[0],len(self.operators['sp']),len(self.operators['sp'])),dtype=complex128)
            for i,ks in enumerate(kmesh):
                result[i,:,:]=self.pt(ks)
            self.cache['pt_mesh']=result
            return result

    def gf_vca(self,omega=None,ks=[]):
        '''
        Returns the single particle Green's function of the system.
        '''
        ngf,ngf_vca,gf=len(self.operators['sp']),len(self.operators['csp']),self.gf(omega)
        return gf_contract(ks=ks,gf_buff=dot(gf,inv(identity(ngf,dtype=complex128)-dot(self.pt(ks),gf))),seqs=self.clmap['seqs'],coords=self.clmap['coords'])/(ngf/ngf_vca)

    def gf_vca_kmesh(self,omega,kmesh):
        '''
        Returns the mesh of the single particle Green's functions of the system.
        '''
        ngf,ngf_vca=len(self.operators['sp']),len(self.operators['csp'])
        gf=self.gf(omega)
        buff=einsum('jk,ikl->ijl',gf,inv(identity(ngf,dtype=complex128)-dot(self.pt_mesh(kmesh),gf)))
        result=zeros((kmesh.shape[0],ngf_vca,ngf_vca),dtype=complex128)
        for n,ks in enumerate(kmesh):
            result[n,:,:]=gf_contract(ks=ks,gf_buff=buff[n,:,:],seqs=self.clmap['seqs'],coords=self.clmap['coords'])
        return result/(ngf/ngf_vca)

def has_integer_solution(coords,vectors):
    nvectors=len(vectors)
    ndim=len(vectors[0])
    a=zeros((3,3))
    for i in xrange(nvectors):
        a[0:ndim,i]=vectors[i]
    if nvectors==2:
        if ndim==2:
            buff=zeros(3)
            buff[2]=cross(vectors[0],vectors[1])
        else:
            buff=cross(vectors[0],vectors[1])
        a[:,2]=buff
    if nvectors==1:
        buff1=a[:,0]
        for i in xrange(3):
            buff2=zeros(3)
            buff2[i]=pi
            if not is_parallel(buff1,buff2): break
        buff3=cross(buff1,buff2)
        a[:,1]=buff2
        a[:,2]=buff3
    b=zeros(3)
    b[0:len(coords)]=coords
    x=dot(inv(a),b)
    if max(abs(x-around(x)))<RZERO:
        return True
    else:
        return False

def VCAEB(engine,app):
    engine.cache.pop('pt_mesh',None)
    erange=linspace(app.emin,app.emax,app.ne)
    result=zeros((app.path.rank,app.ne))
    for i,omega in enumerate(erange):
        result[:,i]=-2*imag((trace(engine.gf_vca_kmesh(omega+engine.mu+app.delta*1j,app.path.mesh),axis1=1,axis2=2)))
    if app.save_data:
        buff=zeros((app.path.rank*app.ne,3))
        for k in xrange(buff.shape[0]):
            i,j=divmod(k,app.path.rank)
            buff[k,0]=j
            buff[k,1]=erange[i]
            buff[k,2]=result[j,i]
        savetxt(engine.dout+'/'+engine.name.full_name+'_EB.dat',buff)
    if app.plot:
        krange=array(xrange(app.path.rank))
        plt.title(engine.name.full_name+'_EB')
        plt.colorbar(plt.pcolormesh(tensordot(krange,ones(app.ne),axes=0),tensordot(ones(app.path.rank),erange,axes=0),result))
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full_name+'_EB.png')
        plt.close()

def VCAFS(engine,app):
    pass

def VCADOS(engine,app):
    engine.cache.pop('pt_mesh',None)
    erange=linspace(app.emin,app.emax,app.ne)
    result=zeros((app.ne,2))
    for i,omega in enumerate(erange):
        result[i,0]=omega
        result[i,1]=-2*imag(sum((trace(engine.gf_vca_kmesh(omega+engine.mu+app.delta*1j,app.BZ.mesh),axis1=1,axis2=2))))
    if app.save_data:
        savetxt(engine.dout+'/'+engine.name.full_name+'_DOS.dat',result)
    if app.plot:
        plt.title(engine.name.full_name+'_DOS')
        plt.plot(result[:,0],result[:,1])
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full_name+'_DOS.png')
        plt.close()
