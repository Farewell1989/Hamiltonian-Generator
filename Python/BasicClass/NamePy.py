class Name(dict):
    '''
    The Name class provides a physical model with a name, whose attributes are as follows:
    1) prefix: a string to describe the physical system, e.g., the lattice in which it is;
    2) suffix: additional remarks on the physical model, e.g., the methods used to calculate physical quantities;
    3) full_name: since the class itself is a list, the parameters of the physical model are included, then the full_name is of the form 'prefix'_'parameters'_'suffix'.
    '''
    
    def __init__(self,prefix='',suffix=''):
        self.prefix=prefix
        self.suffix=suffix
    
    def __str__(self):
        return self.full_name

    @property
    def full_name(self):
        result=self.prefix+'_'
        for obj in self.itervalues():
            result+=repr(obj)+'_'
        result+=self.suffix
        return result

# The following codes are used for tests only.
def test_name():
    a=Name('Hexagon','CPT')
    a[0]=1.0
    a[1]=2.0+2.0j
    print a
    del a[0]
    print a