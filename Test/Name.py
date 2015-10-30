from BasicClass.NamePy import *
def test_name():
    a=Name('Hexagon','CPT')
    a[0]=1.0
    a[1]=2.0+2.0j
    print a
    del a[0]
    print a
