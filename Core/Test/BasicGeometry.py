from Hamiltonian.Core.BasicClass.BasicGeometryPy import *
def test_basicgeometry():
    test_basicgeometry_functions()
    test_point()

def test_basicgeometry_functions():
    a=array([1.0,-1.0,0.0])
    print azimuthd(a)
    print azimuth(a)
    b=array([1.0,1.0,0.0])
    print inner(a,b)
    print cross(a,b)
    c=array([1.0,0.0,0.0])
    print volume(a,b,c)
    print polar(array([1.0,0.0,0.0]))
    print polard(array([1.0,0.0,-1.0]))
    print is_parallel(a,b)

def test_point():
    a=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=Fermi(norbital=2))
    b=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=Fermi(norbital=2,nspin=2,nnambu=1))
    print b
    print 'b==a: ',b==a
    print 'b is a: ',b is a
    table=b.table(priority=lambda key: key.to_tuple(indication='PNSCO'))
    print table
    for i in xrange(len(table)):
        c=b.struct.state_index(i)
        print 'seq_state: ',b.struct.seq_state(**c)
        print c
