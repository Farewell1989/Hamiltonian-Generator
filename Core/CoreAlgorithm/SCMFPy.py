'''
Simple self-consistent mean field theory.
'''
from TBAPy import *
class SCMF(TBA):
    '''
    '''
    def __init__(self,filling=0,mu=0,lattice=None,terms=None,orders=None,nambu=False,**karg):
        '''
        Constructor.
        '''
        self.filling=filling
        self.mu=mu
        self.lattice=lattice
        self.terms=terms
        self.generators={}
        self.generators['h']=Generator(bonds=lattice.bonds,table=lattice.table(nambu),terms=terms,nambu=nambu,half=True)
        self.generators['op']=Generator(bonds=lattice.bonds,table=lattice.table(nambu),terms=orders,nambu=nambu,half=True)
        self.name.update(const=self.generators['h'].parameters['const'])
        self.name.update(alter=self.generators['h'].parameters['alter'])
        #self.name.update(const=self.generators['op'].parameters['const'])
        self.name.update(alter=self.generators['op'].parameters['alter'])
        self.op_m={}
        self.set_op_m()

    def set_op_m(self):
        nmatrix=len(self.generators['h'].table)
        for key in self.generators['op'].terms['alter'].keys():
            self.op_m[key]=zeros((nmatrix,nmatrix),dtype=complex128)
            for opt in self.generators['op'].operators:
                self.op_m[key][opt.seqs]+=opt.value
            self.op_m[key]+=conjugate(self.op_m[key].T)
