from TermPy import *
from IndexPackagePy import *
from IndexPy import *
from BondPy import *
from TablePy import *
from OperatorPy import * 
class Quadratic(Term):
    '''
    Class Quadratic provides a complete and unified description for quadratic terms, i.e. hopping terms, onsite terms and pairing terms. It has the following attributes:
    1) mode: a string to distinguish which kind of quadratic terms it is;
    2) value: the overall coefficient for a quadratic term;
    3) neighbour: as its literal meaning;
    4) indexpackages: two cases,
        (1) an instance of IndexPackageList, which can be referenced in the module IndexPackagePy;
        (2) a function designed to handle those cases where the indexpackages are non-isotropic.
    5) amplitude : a function designed to deal with non-isotropic cases with a real or complex returned value, which depends on the bond on which the quadratic term is defined, and will be multiplied in addition to the overall coefficient value.
    '''
    
    def __init__(self,mode,tag,value,neighbour=0,atoms=[],orbitals=[],spins=[],indexpackages=None,amplitude=None,modulate=None):
        '''
        Note: atoms,orbitals and spins will be packed by an instance of IndexPackageList.
        '''
        super(Quadratic,self).__init__(mode,tag,value,modulate)
        self.neighbour=neighbour
        if not indexpackages is None:
            if isinstance(indexpackages,IndexPackageList):
                self.indexpackages=IndexPackage(1,atoms=atoms,orbitals=orbitals,spins=spins)*indexpackages
            elif callable(indexpackages):
                self.indexpackages=indexpackages
            else:
                raise ValueError('Quadratic init error: the input indexpackages should be an instance of IndexPackageList or a function.')
        else:
            self.indexpackages=IndexPackageList(IndexPackage(1,atoms=atoms,orbitals=orbitals,spins=spins))
        self.amplitude=amplitude

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result='Mode & tag: '+self.mode+' & '+self.tag+'\n'+'Value: '+str(self.value)+'\n'
        if isinstance(self.indexpackages,IndexPackageList):
            result+=str(self.indexpackages)
        if hasattr(self,'modulate'):
            result+='Modulate function: '+str(self.modulate)+'\n'
        return result

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a Quadratic instance with a Quadratic/QuadraticList instance.
        '''
        result=QuadraticList(self.mode)
        result.append(deepcopy(self))
        if isinstance(other,Quadratic):
            result.append(deepcopy(other))
        elif isinstance(other,QuadraticList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('Quadratic "+" error: the other parameter must be an instance of Quadratic or QuadraticList.')
        return result

    def __pos__(self):
        '''
        Overloaded operator(+), i.e. +self.
        '''
        result=QuadraticList()
        result.append(deepcopy(self))
        return result

    def mesh(self,bond,half=True,dtype=complex128):
        '''
        Generate the mesh of a quadratic term defined on a bond.
        '''
        n1=bond.epoint.norbital*bond.epoint.nspin*bond.epoint.nnambu
        n2=bond.spoint.norbital*bond.spoint.nspin*bond.spoint.nnambu
        result=zeros((n1,n2),dtype=dtype)
        if self.neighbour==bond.neighbour:
            value=self.value*(1 if self.amplitude==None else self.amplitude(bond))
            if callable(self.indexpackages):
                buff=self.indexpackages(bond)
            else:
                buff=self.indexpackages
            for obj in buff:
                pv=value*obj.value
                eatom=bond.epoint.atom
                satom=bond.spoint.atom
                if hasattr(obj,'atoms'):
                    eatom=obj.atoms[0]
                    satom=obj.atoms[1]
                if eatom==bond.epoint.atom and satom==bond.spoint.atom:
                    enambu=CREATION if self.mode=='pr' else ANNIHILATION
                    snambu=ANNIHILATION
                    if hasattr(obj,'spins'):
                        if hasattr(obj,'orbitals'):
                            i=bond.epoint.seq_state(obj.orbitals[0],obj.spins[0],enambu)
                            j=bond.spoint.seq_state(obj.orbitals[1],obj.spins[1],snambu)
                            result[i,j]+=pv
                        elif bond.epoint.norbital==bond.spoint.norbital:
                            for k in xrange(bond.epoint.norbital):
                                i=bond.epoint.seq_state(k,obj.spins[0],enambu)
                                j=bond.spoint.seq_state(k,obj.spins[1],snambu)
                                result[i,j]+=pv
                        else:
                            raise ValueError("Quadratic mesh error: the norbital of epoint and spoint of the input bond should be equal.")
                    elif bond.epoint.nspin==bond.spoint.nspin:
                        if hasattr(obj,'orbitals'):
                            for k in xrange(bond.epoint.nspin):
                                i=bond.epoint.seq_state(obj.orbitals[0],k,enambu)
                                j=bond.spoint.seq_state(obj.orbitals[1],k,snambu)
                                result[i,j]+=pv
                        elif n1==n2:
                            if self.mode=='pr':
                                ns=bond.epoint.norbital*bond.epoint.nspin
                                for k in xrange(ns):
                                    result[k,k+ns]+=pv
                            else:
                                for k in xrange(bond.epoint.norbital*bond.epoint.nspin):
                                    result[k,k]+=pv
                        else:
                            raise ValueError("Quadratic mesh error: the norbital of epoint and spoint of the input bond should be equal.")
                    else:
                        raise ValueError("Quadratic mesh error: the nspin of epoint and spoint of the input bond should be equal.")
            if self.neighbour==0:
                for i in xrange(n1):
                    if half: result[i,i]/=2
                    for j in xrange(i):
                        if abs(result[i,j]-conjugate(result[j,i]))<RZERO: result[i,j]=0
        return result

def Hopping(tag,value,neighbour=1,atoms=[],orbitals=[],spins=[],indexpackages=None,amplitude=None,modulate=None):
    '''
    A specified function to construct a hopping term.
    '''
    return Quadratic('hp',tag,value,neighbour,atoms,orbitals,spins,indexpackages,amplitude,modulate)

def Onsite(tag,value,neighbour=0,atoms=[],orbitals=[],spins=[],indexpackages=None,amplitude=None,modulate=None):
    '''
    A specified function to construct an onsite term.
    '''
    return Quadratic('st',tag,value,neighbour,atoms,orbitals,spins,indexpackages,amplitude,modulate)

def Pairing(tag,value,neighbour=0,atoms=[],orbitals=[],spins=[],indexpackages=None,amplitude=None,modulate=None):
    '''
    A specified function to construct an pairing term.
    '''
    return Quadratic('pr',tag,value,neighbour,atoms,orbitals,spins,indexpackages,amplitude,modulate)

class QuadraticList(list):
    '''
    The QuadraticList class pack several Quadratic instances as a whole for convenience.
    '''
    
    def __init__(self,*arg):
        for obj in arg:
            if isinstance(obj,Quadratic):
                self.append(obj)
            else:
                raise ValueError("QuadraticList init error: the input argument should be Quadratic instances.")

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result=''
        for i,v in enumerate(self):
            result+='Node '+str(i)+':\n'+str(v)
        return result

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of a QuadraticList instance with a Quadratic/QuadraticList instance.
        '''
        result=deepcopy(self)
        if isinstance(other,Quadratic):
            result.append(deepcopy(other))
        elif isinstance(other,QuadraticList):
            result.extend(deepcopy(other))
        else:
            raise ValueError('QuadraticList "+" error: the other parameter must be an instance of Quadratic or QuadraticList.')
        return result

    def __mul__(self,other):
        '''
        Overloaded operator(*), which supports the left multiplication with a scalar.
        '''
        result=QuadraticList()
        for obj in self:
            result.append(obj*other)
        return result

    def __rmul__(self,other):
        '''
        Overloaded operator(*), which supports the right multiplication with a scalar.
        '''
        return self.__mul__(other)

    def mesh(self,bond,half,mask=None,dtype=complex128):
        '''
        Generate the mesh of all quadratic terms defined on a bond.
        '''
        if bond.epoint.nnambu==bond.spoint.nnambu:
            n1=bond.epoint.norbital*bond.epoint.nspin*bond.epoint.nnambu
            n2=bond.spoint.norbital*bond.spoint.nspin*bond.spoint.nnambu
            result=zeros((n1,n2),dtype=dtype)
            for obj in self:
                if mask is None or mask(obj):
                    result+=obj.mesh(bond,half,dtype=dtype)
            return result
        else:
            raise ValueError('Quadratic mesh error: the nnambu of epoint and spoint must be equal.')

    def operators(self,bond,table,half=True,dtype=complex128):
        '''
        Generate the set of non-zero operators defined on the input bond.
        1) The index sequences are determined by the index-sequence table.
        2) Because of the hermiticity of the Hamiltonian, when the parameter 'half' is set to be true, only one half of the set is returned. Note that the coefficient of the self-hermitian operators is also divided by a factor 2 so that the whole set exactly equals the returned set plus its Hermitian conjugate.
        3) As for the BdG case, when half=True, only the electron part of the hopping terms and onsite terms are generated since the hole part is nothing but the minus electron part in the matrix representation. As a result, only one quarter of the hopping terms and onsite terms are returned.
        '''
        result=_operators(self.mesh(bond,half,dtype=dtype),bond,table,half)
        if bond.neighbour!=0:
            result.extend(_operators(self.mesh(bond.reversed(),half,mask=lambda quadratic: True if quadratic.mode=='pr' else False,dtype=dtype),bond.reversed(),table,half))
        return result

def _operators(mesh,bond,table,half=True):
    result=OperatorList()
    indices=argwhere(abs(mesh)>RZERO)
    for (i,j) in indices:
        eindex=Index(scope=bond.epoint.scope,site=bond.epoint.site,**bond.epoint.state_index(i))
        sindex=Index(scope=bond.spoint.scope,site=bond.spoint.site,**bond.epoint.state_index(j))
        if eindex in table and sindex in table:
            result.append(E_Quadratic(mesh[i,j],indices=deepcopy([eindex.dagger,sindex]),rcoords=[bond.rcoord],icoords=[bond.icoord],seqs=[table[eindex],table[sindex]]))
            if not half and eindex!=sindex:
                result.append(E_Quadratic(conjugate(mesh[i,j]),indices=deepcopy([sindex.dagger,eindex]),rcoords=[-bond.rcoord],icoords=[-bond.icoord],seqs=[table[sindex],table[eindex]]))
    return result
