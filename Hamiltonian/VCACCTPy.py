from VCAPy import *
class VCACCT(VCA):
    '''
    '''
    def __init__(self,name=None,ensemble='c',filling=0.5,mu=0,nspin=1,cell=None,lattice=None,subsystems=None,terms=None,weiss=None,nambu=False,**karg):
        self.name=Name(prefix=name,suffix=self.__class__.__name__)
        self.ensemble=ensemble
        self.filling=filling
        self.mu=mu
        if self.ensemble.lower()=='c':
            self.name['filling']=self.filling
        elif self.ensemble.lower()=='g':
            self.name['mu']=self.mu
        self.cell=cell
        self.lattice=lattice
        self.terms=terms
        self.weiss=weiss
        self.nambu=nambu
        self.generators={}
        self.generators['pt']=Generator(
                    bonds=      [bond for bond in lattice.bonds if not bond.is_intra_cell() or bond.spoint.scope!=bond.epoint.scope],
                    table=      Table(lattice.indices(nambu=nambu)),
                    terms=      terms if weiss is None else terms+[term*(-1) for term in weiss],
                    nambu=      nambu,
                    half=       True
                    )
        self.name.update(self.generators['pt'].parameters['const'])
        self.name.update(self.generators['pt'].parameters['alter'])
        self.subsystems={}
        for i,(basis,lattice) in enumerate(subsystems):
            self.subsystems[lattice.name]=ONR(
                    name=       lattice.name,
                    ensemble=   ensemble,
                    filling=    filling,
                    mu=         mu,
                    basis=      basis,
                    nspin=      nspin,
                    lattice=    lattice,
                    terms=      terms if weiss is None else terms+weiss,
                    nambu=      nambu,
                    **karg
                )
            if i==0: flag=self.subsystems[lattice.name].nspin
            if flag!=self.subsystems[lattice.name].nspin:
                raise ValueError("VCACCT init error: all the subsystems must have the same nspin.")
        self.nspin=flag
        self.operators={}
        self.set_operators()
        self.clmap={}
        self.set_clmap()
        self.cache={}

    def set_operators(self):
        self.set_operators_perturbation()
        self.set_operators_single_particle()
        self.set_operators_cell_single_particle()

    def gf(self,omega=None):
        pass

def VCACCTGFC(engine,app):
    pass

def VCACCTGF(engine,app):
    pass












