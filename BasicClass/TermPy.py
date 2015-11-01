from numpy import *
from copy import deepcopy
class Term(object):
    '''
    Class Term is the base class for all kinds of terms contained in a Hamiltonian. It has the following attributes:
    1) mode: the type of the term;
    2) tag: the tag specifying the term used for dictionary lookup;
    3) value: the overall coefficient(s) of the term;
    4) modulate: a function used to alter the value of the term.
    '''
    def __init__(self,mode,tag,value,modulate=None):
        self.mode=mode
        self.tag=tag
        self.value=array(value)
        if not modulate is None:
            self.modulate=modulate

    def __mul__(self,other):
        '''
        Overloaded operator(*) which supports the left multiplication with a scalar.
        '''
        result=deepcopy(self)
        result.value*=other
        return result
 
    def __rmul__(self,other):
        '''
        Overloaded operator(*) which supports the right multiplication with a scalar.
        '''
        return self.__mul__(other)
