from VCAPy import *
class VCACCT(VCA):
    '''
    '''
    def __init__(self,name=None,ensemble='c',filling=0.5,mu=0,nspin=1,cell=None,lattice=None,subsystems=None,terms=None,weiss=None,nambu=False,**karg):
        self.name=Name(prefix=name,suffix=self.__class__.__name__)
        self.ensemble=ensemble
        self.filling=filling
        self.mu=mu
        self.cell=cell
        if self.ensemble.lower()=='c':
            self.name['filling']=self.filling
        elif self.ensemble.lower()=='g':
            self.name['mu']=self.mu
        self.generator=Generator(lattice=lattice,terms=terms,nambu=nambu)
        self.name.update(self.generator.parameters['const'])
        self.name.update(self.generator.parameters['alter'])
        self.subsystems={}
        for i,(basis,lattice) in enumerate(subsystems):
            self.subsystems[lattice.name]=ONR(
                    name=       lattice.name,
                    ensemble=   ensemble,
                    filling=    filling,
                    mu=         mu,
                    basis=      basis,
                    nspin=      nspin,
                    generator=  Generator(
                        lattice=    lattice,
                        terms=      terms if weiss is None else terms+weiss,
                        nambu=      nambu
                        ),
                    **karg
                )
            if i==0: flag=self.subsystems[lattice.name].nspin
            if flag!=self.subsystems[lattice.name].nspin:
                raise ValueError("VCACCT init error: all the subsystems must have the same nspin.")
        self.nspin=flag
        self.weiss=Generator(lattice=lattice,terms=weiss,nambu=nambu)
        self.name.update(self.weiss.parameters['const'])
        self.name.update(self.weiss.parameters['alter'])
        self.operators={}
        #self.set_operators()
        self.clmap={}
        #self.set_clmap()
        self.cache={}

    def set_operators(self):
        self.set_operators_perturbation()
        self.set_operators_single_particle()
        self.set_oeprators_cell_single_particle()

    def set_operators_perturbation(self):
        self.operators['pt']=OperatorList()
        buff=self.generator.operators
        table=self.generator.table if self.nspin==2 else subset(self.generator.table,mask=lambda index: True if index.spin==0 else False)
        for opt in buff:
            if norm(opt.icoords)>RZERO:
                if opt.indices[1] in table:
                    self.operators['pt'].append(opt)
        if not self.weiss is None:
            buff=self.weiss.operators
            for opt in buff:
                if norm(opt.icoords)>RZERO:
                    if opt.indices[1] in table:
                        self.operators['pt'].append(opt*(-1))

    def gf(self,omega=None):
        pass


def VCACCTGFC(engine,app):
    pass

def VCACCTGF(engine,app):
    pass












