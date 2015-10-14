from ONRPy import *
from BasicClass.BaseSpacePy import *
class CPT(ONR):
    '''
    '''
    def __init__(self,name=None,ensemble='c',filling=0.5,mu=0,basis=None,nspin=1,cell=None,generator=None,**karg):
        self.cell=cell
        self.ctable=Table(self.cell.indices(nambu=generator.nambu))
        super(CPT,self).__init__(
            name=       name,
            ensemble=   ensemble,
            filling=    filling,
            mu=         mu,
            basis=      basis,
            nspin=      nspin,
            generator=  generator
            )
        self.clmap=[]
        self.set_clmap()
        self.cache={}

    def set_operators(self):
        self.set_operators_hamiltonian_and_perturbation()
        self.set_operators_single_particle()
        self.set_operators_cell_single_particle()

    def set_operators_hamiltonian_and_perturbation(self):
        self.operators['h']=OperatorList()
        self.operators['pt']=OperatorList()
        buff=self.generator.operators
        table=self.generator.table if self.nspin==2 else subset(self.generator.table,mask=lambda index: True if index.spin==0 else False)
        for opt in buff:
            if norm(opt.icoords)>RZERO:
                if opt.indices[1] in table:
                    self.operators['pt'].append(opt)
            else:
                self.operators['h'].append(opt)

    def set_operators_cell_single_particle(self):
        self.operators['csp']=OperatorList()
        table=self.ctable if self.nspin==2 else subset(self.ctable,mask=lambda index: True if index.spin==0 else False)
        for index,seq in table.iteritems():
            if isinstance(index,Index): 
                self.operators['csp'].append(E_Linear(1,indices=[index],rcoords=[self.cell.points[index.site].rcoord],icoords=[self.cell.points[index.site].icoord],seqs=[seq]))
        self.operators['csp'].sort(key=lambda opt: opt.seqs[0])

    def set_clmap(self):
        for i in xrange(len(self.operators['csp'])):
            self.clmap.append(OperatorList())
        for optl in self.operators['sp']:
            for i,optc in enumerate(self.operators['csp']):
                if optl.indices[0].orbital==optc.indices[0].orbital and optl.indices[0].spin==optc.indices[0].spin and optl.indices[0].nambu==optc.indices[0].nambu and has_integer_solution(optl.rcoords[0]-optc.rcoords[0],self.cell.vectors):
                    self.clmap[i].append(optl)
                    break

    def gf(self,omega):
        if not 'GF' in self.apps:
            self.addapps(app=GF((len(self.operators['sp']),len(self.operators['sp'])),run=ONRGF))
        self.apps['GF'].omega=omega
        self.runapps('GF')
        return self.apps['GF'].gf

    def gf_mesh(self,omegas):
        if 'gf_mesh' in self.cache:
            return self.cache['gf_mesh']
        else:
            result=zeros((omegas.shape[0],len(self.operators['sp']),len(self.operators['sp'])),dtype=complex128)
            for i,omega in enumerate(omegas):
                result[i,:,:]=self.gf(omega)
            self.cache['gf_mesh']=result
            return result

    def pt(self,ks):
        ngf=len(self.operators['sp'])
        result=zeros((ngf,ngf),dtype=complex128)
        for opt in self.operators['pt']:
            result[opt.seqs]+=opt.value*(1 if len(ks)==0 else exp(-1j*inner(ks,opt.icoords[0])))
        return result+conjugate(result.T)

    def pt_mesh(self,kmesh):
        if 'pt_mesh' in self.cache:
            return self.cache['pt_mesh']
        else:
            result=zeros((kmesh.shape[0],len(self.operators['sp']),len(self.operators['sp'])),dtype=complex128)
            for i,ks in enumerate(kmesh):
                result[i,:,:]=self.pt(ks)
            self.cache['pt_mesh']=result
            return result

    def gf_cpt(self,omega=None,ks=None):
        ngf=len(self.operators['sp'])
        ngf_cpt=len(self.operators['csp'])
        result=zeros((ngf_cpt,ngf_cpt),dtype=complex128)
        gf=self.apps['GF'].gf if omega==None else self.gf(omega)
        buff=dot(gf,inv(identity(ngf,dtype=complex128)-dot(self.pt(ks),gf)))
        for i in xrange(ngf_cpt):
            for k,optk in enumerate(self.clmap[i]):
                row=optk.seqs[0]
                for j in xrange(ngf_cpt):
                    for l,optl in enumerate(self.clmap[j]):
                        col=optl.seqs[0]
                        result[i,j]+=buff[row,col]*exp(1j*inner(ks,optl.rcoords[0]-optk.rcoords[0]))
        return result/(ngf/ngf_cpt)

#    def cl_mesh(self,kmesh):
#        if 'cl_mesh' in self.cache:
#            return self.cache['cl_mesh']
#        else:
#            nmesh=max([len(opts) for opts in self.clmap])
#            result=zeros((kmesh.shape[0],nmesh,nmesh),dtype=complex128)
            
    def gf_cpt_kmesh(self,omega,kmesh):
        import pdb
        ngf=len(self.operators['sp'])
        gf=self.gf(omega)
        buff=einsum('jk,ikl->ijl',gf,inv(identity(ngf,dtype=complex128)-dot(self.pt_mesh(kmesh),gf)))
        ngf_cpt=len(self.operators['csp'])
        result=zeros((kmesh.shape[0],ngf_cpt,ngf_cpt),dtype=complex128)
        for n,ks in enumerate(kmesh):
            for i in xrange(ngf_cpt):
                for k,optk in enumerate(self.clmap[i]):
                    row=optk.seqs[0]
                    for j in xrange(ngf_cpt):
                        for l,optl in enumerate(self.clmap[j]):
                            col=optl.seqs[0]
                            result[n,i,j]+=buff[n,row,col]*exp(1j*inner(ks,optl.rcoords[0]-optk.rcoords[0]))
        return result/(ngf/ngf_cpt)

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

def CPTEB(engine,app):
    erange=linspace(app.emin,app.emax,app.ne)
    result=zeros((app.path.rank,app.ne))
    engine.cache.pop('pt_mesh',None)
    for i,omega in enumerate(erange):
        result[:,i]=-2*imag((trace(engine.gf_cpt_kmesh(omega+engine.mu+app.delta*1j,app.path.mesh),axis1=1,axis2=2)))
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
        plt.pcolor(tensordot(krange,ones(app.ne),axes=0),tensordot(ones(app.path.rank),erange,axes=0),result)
        if app.show:
            plt.show()
        else:
            plt.savefig(engine.dout+'/'+engine.name.full_name+'_EB.png')
        plt.close()

def CPTFS(engine,app):
    pass

def CPTDOS(engine,app):
    engine.addapps(app=GF((len(engine.operators['sp']),len(engine.operators['sp'])),run=ONRGF))
    erange=linspace(app.emin,app.emax,app.ne)
    result=zeros((app.ne,2))
    engine.cache.pop('pt_mesh',None)
    for i,omega in enumerate(erange):
        result[i,0]=omega
        result[i,1]=-2*imag(sum((trace(engine.gf_cpt_kmesh(omega+engine.mu+app.delta*1j,app.BZ.mesh),axis1=1,axis2=2))))
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

# The following codes are used for tests only.
def test_cpt():
    U=0.0
    t=-1.0
    m=2;n=2
    p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],atom=0,norbital=1,nspin=2,nnambu=1,scope='WG'+str(m)+str(n))
    a1=array([1.0,0.0])
    a2=array([0.0,1.0])
    a=CPT(
            din=        'Results/Coeff',
            dout=       'Results',
            name=       'WG'+str(m)+str(n),
            ensemble=   'c',
            filling=    0.5,
            mu=         U/2,
            basis=      BasisE(up=(m*n,m*n/2),down=(m*n,m*n/2)),
            #basis=      BasisE((2*m*n,m*n)),
            nspin=      1,
            cell=       Lattice(name='WG',points=[p1],vectors=[a1,a2]),
            generator=  Generator(
                lattice=    Lattice(name='WG'+str(m)+str(n),points=[p1],translations=[(a1,m),(a2,n)],vectors=[a1*m,a2*n]),
                terms=[     Hopping('t',t,neighbour=1),
                            Hubbard('U',[U])
                            ]
                )
            )
    a.addapps('GFC',GFC(nstep=200,save_data=False,vtype='SY',run=ONRGFC))
    a.addapps('DOS',DOS(BZ=square_bz(nk=50),emin=-5,emax=5,ne=400,delta=0.05,save_data=False,run=CPTDOS,plot=True,show=False))
    a.addapps('EB',EB(path=square_gxm(nk=100),emax=6.0,emin=-6.0,delta=0.05,ne=400,save_data=False,plot=True,show=True,run=CPTEB))
    a.runapps()
