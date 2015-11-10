from collections import OrderedDict
class Name:
    '''
    The Name class provides a physical model with a name, whose attributes are as follows:
    1) prefix: a string to describe the physical system, e.g., the lattice in which it is;
    2) suffix: additional remarks on the physical model, e.g., the methods used to calculate physical quantities;
    3) _alter: an ordered dict containing the contant parameters;
    4) _const: an ordered dict containing the alterable parameters.
    '''
    
    def __init__(self,prefix='',suffix=''):
        self.prefix=prefix
        self.suffix=suffix
        self._const=OrderedDict()
        self._alter=OrderedDict()
    
    def __str__(self):
        return self.full

    def update(self,const=None,alter=None):
        if const is not None:
            self._const.update(const)
        if alter is not None:
            self._alter.update(alter)

    @property
    def const(self):
        '''
        Return the name string containing only contant parameters.
        '''
        result=self.prefix+'_'
        for obj in self._const.itervalues():
            result+=repr(obj)+'_'
        result+=self.suffix
        return result

    @property
    def alter(self):
        '''
        Return the name string containing only alterable parameters.
        '''
        result=self.prefix+'_'
        for obj in self._alter.itervalues():
            result+=repr(obj)+'_'
        result+=self.suffix
        return result

    @property
    def full(self):
        '''
        Return the name string containing both contant parameters and alterable parameters.
        '''
        result=self.prefix+'_'
        for obj in self._const.itervalues():
            result+=repr(obj)+'_'
        for obj in self._alter.itervalues():
            result+=repr(obj)+'_'
        result+=self.suffix
        return result
