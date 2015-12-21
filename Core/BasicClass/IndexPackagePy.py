'''
IndexPackage.
'''
from numpy import *
from copy import deepcopy

class IndexPackage:
    '''
    This class assumes part of a systematic description of a general quadratic term.
    Attributes:
        value: float or complex
            The overall coefficient of the index pack.
        atoms: 1D ndarray of integers with len==2
            The atom indices for the quadratic term.
        orbitals: 1D ndarray of integers with len==2
            The orbital indices for the quadratic term.
        spins: 1D ndarray of integers with len==2
            The spin indices for the quadratic term.
    '''
    
    def __init__(self,value,atom1=None,atom2=None,orbital1=None,orbital2=None,spin1=None,spin2=None,atoms=None,orbitals=None,spins=None):
        '''
        Constructor. 
        It can be used in two different ways:
        1) IndexPackage(value,atom1=...,atom2=...,orbital1=...,orbital2=...,spin1=...,spin2=...)
        2) IndexPackage(value,atoms=...,orbitals=...,spins=...)
        Parameters:
            value: float or complex
                The overall coefficient of the index pack
            atom1,atom2: integer,optional
                The atom indices.
            orbital1,orbital2: integer,optional
                The orbital indices.
            spin1,spin2: integer, optional
                The spin indices.
            atoms: 1D array-like of integers with len==1,2,optional
                The atom indices.
            orbitals: 1D array-like of integers with len==1,2,optional
                The orbital indices.
            spins: 1D array-like of integers with len==1,2,optional
                The spin indices.
        '''
        self.value=value
        if atom1 is not None and atom2 is not None: self.atoms=array([atom1,atom2])
        if orbital1 is not None and orbital2 is not None: self.orbitals=array([orbital1,orbital2])
        if spin1 is not None and spin2 is not None: self.spins=array([spin1,spin2])
        if atoms is not None:
            if len(atoms)==2: self.atoms=array(atoms)
            elif len(atoms)==1: self.atoms=array([atoms[0],atoms[0]])
        if orbitals is not None:
            if len(orbitals)==2: self.orbitals=array(orbitals)
            elif len(orbitals)==1: self.orbitals=array([orbitals[0],orbitals[0]])
        if spins is not None:
            if len(spins)==2: self.spins=array(spins)
            elif len(spins)==1: self.spins=array([spins[0],spins[0]])

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
        Overloaded operator(*), which supports the multiplication of an IndexPackage instance with an IndexPackage/IndexPackageList instance or a scalar.
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
    This class packs several IndexPackage as a whole for convenience.
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
        Overloaded operator(*), which supports the multiplication of an IndexPackageList instance with an IndexPackage/IndexPackageList instance or a scalar.
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
