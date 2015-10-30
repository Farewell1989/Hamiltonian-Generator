from numpy import *
from copy import deepcopy

class IndexPackage:
    '''
    The IndexPackage class assumes part of a systematic description of a general quadratic term, e.g. a hopping term, an onsite term or a pairing term. It has the following attributes:
    1) value: the overall coefficient of the index pack;
    2) atoms: a length=2 array containing the atom indices for the quadratic term;
    3) orbitals: a length=2 array containing the orbital indices for the quadratic term;
    4) spins: a length=2 array containing the spin indices for the quadratic term.
    '''
    
    def __init__(self,value,atom1=None,atom2=None,orbital1=None,orbital2=None,spin1=None,spin2=None,atoms=[],orbitals=[],spins=[]):
        self.value=value
        if atom1!=None and atom2!=None: self.atoms=array([atom1,atom2])
        if orbital1!=None and orbital2!=None: self.orbitals=array([orbital1,orbital2])
        if spin1!=None and spin2!=None: self.spins=array([spin1,spin2])
        if len(atoms)==2: self.atoms=array(atoms)
        if len(atoms)==1: self.atoms=array([atoms,atoms])
        if len(orbitals)==2: self.orbitals=array(orbitals)
        if len(orbitals)==1: self.orbitals=array([orbitals,orbitals])
        if len(spins)==2: self.spins=array(spins)
        if len(spins)==1: self.spins=array([spins,spins])

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result='Value:'+str(self.value)+'\n'
        if hasattr(self,'atoms'):result+='Atoms: '+str(self.atoms)+'\n'
        if hasattr(self,'orbitals'): result+='Orbitals: '+str(self.orbitals)+'\n'
        if hasattr(self,'spins'): result+='Spins:'+str(self.spins)+'\n'
        return result

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of an IndexPackage instance with an IndexPackage/IndexPackageList instance.
        '''
        result=IndexPackageList()
        result.append(deepcopy(self))
        if isinstance(other,IndexPackage):
            result.append(deepcopy(other))
        elif isinstance(other,IndexPackageList):
            result.extend(deepcopy(other))
        else:
            raise ValueError("IndexPackage '+' error: the 'other' parameter must be of class IndexPackage or IndexPackageList.")
        return result

    def _mul(self,other):
        '''
        Private methods used for operator(*) overloading.
        '''
        if isinstance(other,IndexPackage):
            result=IndexPackage(self.value*other.value)
            if hasattr(self,'atoms'): result.atoms=copy(self.atoms)
            if hasattr(self,'orbitals'): result.orbitals=copy(self.orbitals)
            if hasattr(self,'spins'): result.spins=copy(self.spins)
            if hasattr(other,'atoms'):
                if not hasattr(result,'atoms'):
                    result.atoms=copy(other.atoms)
                else:
                    raise ValueError("IndexPackage '*' error: 'self' and 'other' cannot simultaneously have the 'atoms' attribute.")
            if hasattr(other,'orbitals'):
                if not hasattr(result,'orbitals'):
                    result.orbitals=copy(other.orbitals)
                else:
                    raise ValueError("IndexPackage '*' error: 'self' and 'other' cannot simultaneously have the 'orbitals' attribute.")
            if hasattr(other,'spins'):
                if not hasattr(result,'spins'):
                    result.spins=copy(other.spins)
                else:
                    raise ValueError("IndexPackage '*' error: 'self' and 'other' cannot simultaneously have the 'spins' attribute.")
        else:
            result=deepcopy(self)
            result.value=self.value*other
        return IndexPackageList(result)
    
    def __mul__(self,other):
        '''
        Overloaded operator(*), which supports the multiplication of an IndexPackage instance with an IndexPackage/IndexPackageList instance.
        '''
        result=IndexPackageList()
        if isinstance(other,IndexPackage) or type(other)==int or type(other)==long or type(other)==float or type(other)==complex:
            result.extend(self._mul(other))
        elif isinstance(other,IndexPackageList):
            for buff in other:
                result.extend(self._mul(buff))
        else:
            raise ValueError("IndexPackage '*' error: the 'other' parameter must be of class IndexPackage or .")
        return result

class IndexPackageList(list):
    '''
    IndexPackageList, which can pack several IndexPackage as a whole for convenience.
    '''
    def __init__(self,*arg):
        for buff in arg:
            if isinstance(buff,IndexPackage):
                self.append(buff)
            else:
                raise ValueError("IndexPackageList init error: the input parameters must be of class IndexPackage.")

    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result=''
        for i,obj in enumerate(self):
            result+='IndexPackage '+str(i)+':\n'+str(obj)
        return result
                
    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of an IndexPackageList instance with an IndexPackage/IndexPackageList instance.
        '''
        result=IndexPackageList(*deepcopy(self))
        if isinstance(other,IndexPackage):
            result.append(deepcopy(other))
        elif isinstance(other,IndexPackageList):
            result.extend(deepcopy(other))
        else:
            raise ValueError("IndexPackageList '+' error: the 'other' parameter must be of class IndexPackage or IndexPackageList.")
        return result
    
    def __mul__(self,other):
        '''
        Overloaded operator(*), which supports the multiplication of an IndexPackageList instance with an IndexPackage/IndexPackageList instance.
        '''
        result=IndexPackageList()
        if isinstance(other,IndexPackage) or isinstance(other,IndexPackageList) or type(other)==int or type(other)==long or type(other)==float or type(other)==complex:
            for buff in self:
                result.extend(buff*other)
        else:
            raise ValueError("IndexPackageList '*' error: the 'other' parameter must be of class IndexPackage or IndexPackageList.")
        return result

def sigmax(mode):
    '''
    The Pauli matrix SigmaX, which can act on the space of spins('sp'), orbitals('ob') or sublattices('sl').
    '''
    result=IndexPackageList()
    if mode.lower()=='sp':
        result.append(IndexPackage(1.0,spin1=0,spin2=1))
        result.append(IndexPackage(1.0,spin1=1,spin2=0))
    elif mode.lower()=='ob':
        result.append(IndexPackage(1.0,orbital1=0,orbital2=1))
        result.append(IndexPackage(1.0,orbital1=1,orbital2=0))
    elif mode.lower()=='sl':
        result.append(IndexPackage(1.0,atom1=0,atom2=1))
        result.append(IndexPackage(1.0,atom1=1,atom2=0))
    else:
        raise ValueError("SigmaX error: mode '"+mode+"' not supported, which must be 'sp', 'ob', or 'sl'.")
    return result

def sigmay(mode):
    '''
    The Pauli matrix SigmaY, which can act on the space of spins('sp'), orbitals('ob') or sublattices('sl').
    '''
    result=IndexPackageList()
    if mode.lower()=='sp':
        result.append(IndexPackage(1.0j,spin1=0,spin2=1))
        result.append(IndexPackage(-1.0j,spin1=1,spin2=0))
    elif mode.lower()=='ob':
        result.append(IndexPackage(1.0j,orbital1=0,orbital2=1))
        result.append(IndexPackage(-1.0j,orbital1=1,orbital2=0))
    elif mode.lower()=='sl':
        result.append(IndexPackage(1.0j,atom1=0,atom2=1))
        result.append(IndexPackage(-1.0j,atom1=1,atom2=0))
    else:
        raise ValueError("SigmaY error: mode '"+mode+"' not supported, which must be 'sp', 'ob', or 'sl'.")
    return result

def sigmaz(mode):
    '''
    The Pauli matrix SigmaZ, which can act on the space of spins('sp'), orbitals('ob') or sublattices('sl').
    '''
    result=IndexPackageList()
    if mode.lower()=='sp':
        result.append(IndexPackage(-1.0,spin1=0,spin2=0))
        result.append(IndexPackage(1.0,spin1=1,spin2=1))
    elif mode.lower()=='ob':
        result.append(IndexPackage(-1.0,orbital1=0,orbital2=0))
        result.append(IndexPackage(1.0,orbital1=1,orbital2=1))
    elif mode.lower()=='sl':
        result.append(IndexPackage(-1.0,atom1=0,atom2=0))
        result.append(IndexPackage(1.0,atom1=1,atom2=1))
    else:
        raise ValueError("SigmaZ error: mode '"+mode+"' not supported, which must be 'sp', 'ob', or 'sl'.")
    return result
