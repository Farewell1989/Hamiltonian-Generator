from GlobalPy import *
from BasicClass.TablePy import *
from BasicClass.OperatorPy import *
from numpy.linalg import norm
class Generator:
    '''
    Class Generator provides methods to generate and update operators of a Hamiltonian according to the lattice of the model and the descriptions of its terms. It has the following attributes:
    1) lattice: the lattice of the model;
    2) parameters: a dict containing all model parameters, divided into two groups, the constant ones and the alterable ones;
    3) terms: a dict containing all terms contained in the model, divided into two groups, the constant ones and the alterable ones;
    4) nambu: a flag to tell whether the nambu space takes part;
    5) table: the index-sequence table of the system;
    6) half: a flag to tell whether the generated operators contains its hermitian conjugate part;
    7) boundary: a function to specify the boundary conditions to control the generation of operators;
    8) cache: the working cache used to handle the generation and update of the operators.
    '''
    def __init__(self,lattice,parameters=None,terms=None,nambu=False,half=True,boundary=None):
        self.lattice=lattice
        self.parameters={}
        self.terms={}
        self.set_parameters_and_terms(parameters,terms)
        self.nambu=nambu
        self.table=Table(self.lattice.indices(nambu=self.nambu))
        self.half=half
        self.boundary=boundary
        self.cache={}
        self.set_cache()

    def set_parameters_and_terms(self,parameters,terms):
        self.parameters['const']={}
        if not parameters is None:
            self.parameters['const'].update(parameters)
        self.parameters['alter']={}
        self.terms['const']={}
        self.terms['alter']={}
        if not terms is None:
            for term in terms:
                if hasattr(term,'modulate'):
                    self.parameters['alter'][term.tag]=term.value
                    self.terms['alter'][term.tag]=+term
                else:
                    self.parameters['const'][term.tag]=term.value
                    lterm=+term
                    if lterm.__class__.__name__ in self.terms['const']:
                        self.terms['const'][lterm.__class__.__name__].extend(lterm)
                    else:
                        self.terms['const'][lterm.__class__.__name__]=lterm

    def set_cache(self):
        if 'const' in self.terms:
            self.cache['const']=OperatorList()
            for bonds in self.lattice.bonds:
                for bond in bonds:
                    if self.boundary is None or self.boundary(bond):
                        for terms in self.terms['const'].itervalues():
                            self.cache['const'].extend(terms.operators(bond,self.table,half=self.half))
        if 'alter' in self.terms:
            self.cache['alter']={}
            for key in self.terms['alter'].iterkeys():
                self.cache['alter'][key]=OperatorList()
                for bonds in self.lattice.bonds:
                    if self.boundary is None or self.boundary(bond):
                        for bond in bonds:
                            self.cache['alter'][key].extend(self.terms['alter'][key].operators(bond,self.table,half=self.half))

    def __str__(self):
        result='Parameters:\n'
        if len(self.parameters['const']>0):
            result+='Constant part: '+str(self.parameters['const'])+'\n'
        if len(self.parameters['alter']>0):
            result+='Alterable part: '+str(self.parameters['alter'])+'\n'
        for key,obj in self.terms['const'].iteritems():
            result+=key+':\n'+str(obj)
        for key,obj in self.terms['alter'].iteritems():
            result+=key+':\n'+str(obj)
        return result

    @property
    def operators(self):
        result=OperatorList()
        if 'const' in self.cache:
            result.extend(self.cache['const'])
        if 'alter' in self.cache:
            for opts in self.cache['alter'].itervalues():
                result.extend(opts)
        return result

    def update(self,**karg):
        if 'alter' in self.terms:
            masks={key:False for key in self.terms['alter'].iterkeys()}
            for key,term in self.terms['alter'].iteritems():
                nv=term[0].modulate(**karg)
                if norm(array(nv)-array(term[0].value))>RZERO:
                    term[0].value=nv
                    self.parameters['alter'][key]=nv
                    masks[key]=True
            for key,mask in masks.iteritems():
                if mask:
                    self.cache['alter'][key]=OperatorList()
                    for bonds in self.lattice.bonds:
                        if self.boundary is None or self.boundary(bond):
                            for bond in bonds:
                                self.cache['alter'][key].extend(self.terms['alter'][key].operators(bond,self.table,half=self.half))