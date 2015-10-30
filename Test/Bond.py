from BasicClass.BondPy import *
def test_bond():
    a=Bond(0,Point(0,rcoord=[0.0,0.0],icoord=[0.0,0.0]),Point(1,rcoord=[0.0,1.0],icoord=[0.0,0.0]))
    print a
    print a.rcoord
    print a.icoord
    print a.is_intra_cell()
