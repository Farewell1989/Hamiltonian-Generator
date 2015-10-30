from BasicClass.BasicGeometryPy import *
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
    a=Point(0,[0.0,0.0],[0.0,0.0])
    b=Point(0,[0.0,0.0],[0.0,0.0],norbital=1,nspin=2,nnambu=1)
    print b
    print 'b==a: ',b==a
    print 'b is a: ',b is a
    for i in xrange(b.norbital*b.nspin*b.nnambu):
        c=b.state_index(i)
        print 'seq_state: ',b.seq_state(**c)
        print c
