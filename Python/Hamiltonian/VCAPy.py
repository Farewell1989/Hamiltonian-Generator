from CPTPy import *
class VCA(CPT):
    '''
    '''
    def __init__(self,name=None,ensemble='c',filling=0.5,mu=0,basis=None,nspin=1,cell=None,generator=None,weiss=None,**karg):
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
        self.cell=cell
        self.ctable=Table(self.cell.indices(nambu=generator.nambu))
        self.weiss=Generator(lattice=self.generator.lattice,terms=weiss,nambu=self.generator.nambu,half=self.generator.half)
        self.name.update(self.weiss.parameters['const'])
        self.name.update(self.weiss.parameters['alter'])
        self.operators={}
        self.set_operators()
        self.clmap=[]
        self.set_clmap()
        self.cache={}

    def set_operators(self):
        self.set_operators_hamiltonian_and_perturbation()
        self.set_operators_single_particle()
        self.set_operators_cell_single_particle()

    def set_operators_hamiltonian_and_perturbation(self):
        super(VCA,self).set_operators_hamiltonian_and_perturbation()
        table=self.generator.table if self.nspin==2 else subset(self.generator.table,mask=lambda index: True if index.spin==0 else False)
        buff=self.weiss.operators
#        print self.weiss
        for opt in buff:
            if norm(opt.icoords)>RZERO:
                if opt.indices[1] in table:
                    self.operators['pt'].append(opt*(-1))
            else:
                self.operators['h'].append(opt)

# The following codes are used for tests only.
def test_vca():
    U=2.0
    t=-1.0
    m=2;n=2
    p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],atom=0,norbital=1,nspin=2,nnambu=1,scope='WG'+str(m)+str(n))
    a1=array([1.0,0.0])
    a2=array([0.0,1.0])
    a=VCA(
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
                ),
            weiss=[     Onsite('afm',0.1,indexpackages=sigmaz('sp'),amplitude=lambda bond: 1 if bond.epoint.site in (0,3) else -1)
                        ]
            )
    a.addapps('GFC',GFC(nstep=200,save_data=False,vtype='SY',run=ONRGFC))
#    a.addapps('DOS',DOS(BZ=square_bz(nk=50),emin=-5,emax=5,ne=400,delta=0.05,save_data=False,run=CPTDOS,plot=True,show=False))
    a.addapps('EB',EB(path=square_gxm(nk=100),emax=6.0,emin=-6.0,delta=0.05,ne=400,save_data=False,plot=True,show=True,run=CPTEB))
    a.runapps()