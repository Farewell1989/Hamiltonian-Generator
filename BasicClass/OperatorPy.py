from IndexPy import *
from BasicGeometryPy import *
from copy import deepcopy
class Operator:
    '''
    The Operator class gives a unified description of different operators with different ranks or types. It has the following attributes:
    1) mode: used to distinguish different types or ranks of operators;
    2) value: the overall coefficient of the operator;
    3) indices: associated indices of the operator, whose length should be equal to the operator's rank;
    4) rcoords: associated real coordinates of the operator;
    5) icoords: associated lattice coordinates of the operator;
    6) seqs: associated sequences of the operator, whose length should be equal to the operator's rank.
    Note: 
    1) The lengths of rcoords and icoords are not forced to be equal to the operator's rank because firstly, some of its rank-1 terms may share the same rcoord or icoord, and secondly, the rcoords and icoords is the whole operator's property instead of each of its rank-1 component. However, for a set of operators with the same attribute mode, the lengths of their rcoords and icoords should be fixed and equal to each other respectively.
    2) Current supported modes include 'e_linear'(rank=1 electron operators), 'e_quadratic'(rank=2 electron operators) and 'e_hubbard'(rank=4 electron operators).For e_quadratic, only one rcoord and icoord is needed because a quadratic operator is always defined on a bond and only the bond's rcoord and icoord is needed. For 'e_hubbard', likewise, only one rcoord and icoord is needed because Hubbard terms are always on-site ones.
    '''
    
    def __init__(self,mode,value,indices,rcoords,icoords,seqs):
        self.mode=mode
        self.value=value
        self.indices=indices
        self.rcoords=[]
        for obj in rcoords:
            self.rcoords.append(array(obj))
        self.icoords=[]
        for obj in icoords:
            self.icoords.append(array(obj))
        self.seqs=tuple(seqs)
    
    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result='Mode: '+self.mode+'\nValue: '+str(self.value)+'\n'
        for i,obj in enumerate(self.indices):
            result+='Index '+str(i)+': '+str(obj)
        for i,obj in enumerate(self.rcoords):
            result+='RCoord '+str(i)+': '+str(obj)+'\n'
        for i,obj in enumerate(self.icoords):
            result+='ICoord '+str(i)+': '+str(obj)+'\n'
        result+='Seqs: '+str(self.seqs)+'\n'
        return result

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of an Operator instance with an Operator/OperatorList instance.
        '''
        result=OperatorList()
        if isinstance(other,Operator):
            if self.is_combinable(other):
                if abs(self.value+other.value)>RZERO:
                    result.append(Opeartor(mode=self.mode,value=self.value+other.value,indices=deepcopy(self.indices),rcoords=self.rcoords,icoords=self.icoords,seqs=self.seqs))
            else:
                result.append(deepcopy(self))
                result.append(deepcopy(other))
        elif isinstance(other,OperatorList):
            mask=True
            for obj in other:
                if mask and self.is_combinable(obj):
                    if abs(self.value+obj.value)>RZERO:
                        result.append(Operator(mode=self.mode,value=self.value+obj.value,indices=deepcopy(self.indices),rcoords=self.rcoords,icoords=self.icoords,seqs=self.seqs))
                    mask=False
                else:
                    result.append(deepcopy(obj))
            if mask: result.append(deepcopy(self))
        else:
            raise ValueError("Operator '+' error: the 'other' parameter must be of class Operator or OperatorList.")
        return result
    
    def __mul__(self,other):
        '''
        Overloaded operator(*), which supports the multiplication of an Operator instance with a scalar.
        '''
        return Operator(mode=self.mode,value=self.value*other,indices=deepcopy(self.indices),rcoords=self.rcoords,icoords=self.icoords,seqs=self.seqs)
        
    def __rmul__(self,other):
        '''
        Overloaded operator(*), which supports the multiplication of an Operator instance with a scalar.
        '''
        return self.__mul__(other)

    def __eq__(self,other):
        '''
        Overloaded operator(==).
        '''
        return abs(self.value-other.value)<RZERO and self.is_combinable(other)
    
    def __ne__(self,other):
        '''
        Overloaded operator(!=).
        '''
        return not self==other
    
    @property
    def dagger(self):
        '''
        The dagger, i.e. the Hermitian conjugate of an operator.
        '''
        indices=[]
        for obj in self.indices:
            indices.append(obj.dagger)
        indices.reverse()
        return Operator(mode=self.mode,value=conjugate(self.value),indices=indices,rcoords=list(reversed(self.rcoords)),icoords=list(reversed(self.icoords)),seqs=fliplr([self.seqs])[0])

    @property
    def rank(self):
        '''
        Return the rank of the operator.
        '''
        return len(self.seqs)

    def is_combinable(self,other):
        '''
        Judge whether two operators are combinable.
        '''
        if self.mode==other.mode and len(self.rcoords)==len(other.rcoords) and len(self.icoords)==len(other.icoords):
            for v1,v2 in zip(self.rcoords,other.rcoords):
                if norm(v1-v2)>RZERO: return False
            for v1,v2 in zip(self.icoords,other.icoords):
                if norm(v1-v2)>RZERO: return False
            if any(self.indices!=other.indices) or any(self.seqs!=other.seqs): return False
            return True
        else:
            return False
        
    def is_normal_ordered(self):
        '''
        Judge whether an operator is normal ordered.
        '''
        buff=True
        for index in self.indices:
            if index.nambu==ANNIHILATION: buff=False
            if not buff and index.nambu==CREATION: return False
        return True

    def is_Hermitian(self):
        '''
        Judge whether an operator is Hermitian.
        '''
        return self==self.dagger

def E_Linear(value,indices,rcoords,icoords,seqs):
    '''
    A specialized constructor to create an Operator instance with mode='e_linear'.
    '''
    return Operator('e_linear',value,indices,rcoords,icoords,seqs)

def E_Quadratic(value,indices,rcoords,icoords,seqs):
    '''
    A specialized constructor to create an Operator instance with mode='e_quadratic'.
    '''
    return Operator('e_quadratic',value,indices,rcoords,icoords,seqs)

def E_Hubbard(value,indices,rcoords,icoords,seqs):
    '''
    A specialized constructor to create an Operator instance with mode='e_hubbard'.
    '''
    return Operator('e_hubbard',value,indices,rcoords,icoords,seqs)

class OperatorList(list):
    '''
    OperatorList, which can pack several operators as a whole for convenience.
    '''
    
    def __init__(self):
        pass
    
    def __str__(self):
        '''
        Convert an instance to string.
        '''
        result=''
        for i,obj in enumerate(self):
           result+='Operator '+str(i)+':\n'+str(obj)
        return result

    def __add__(self,other):
        '''
        Overloaded operator(+), which supports the addition of an OperatorList instance with an Operator/OperatorList instance.
        '''
        if isinstance(other,Operator):
            return other.__add__(self)
        elif isinstance(other,OperatorList):
            result=deepcopy(self)
            for obj in other:
                result=obj.__add__(result)
        else:
            raise ValueError("OperatorList '+' error: the 'other' parameter must be of class Operator or OperatorList.")
        return result

    def __mul__(self,other):
        '''
        Overloaded operator(*), which supports the multiplication of an OperatorList instance with a scalar.
        '''
        result=OperatorList()
        for obj in self:
            result.append(obj*other)
        return result

    def __rmul__(self,other):
        '''
        Overloaded operator(*), which supports the multiplication of an OperatorList instance with a scalar.
        '''
        return self.__mul__(other)
