from BasicClass.LatticePy import *
import time
def test_lattice():
    test_lattice_body() 
    test_lattice_indices()
    #test_super_lattice()   
    
def test_lattice_body():
    p1=Point(site=0,norbital=2,nspin=3,nnambu=2,rcoord=[0.0,0.0],icoord=[0,0])
    a1=array([1.0,0.0])
    a2=array([0.0,1.0])
    m=2;n=2
    stime=time.time()
    a=Lattice('L'+str(m)+str(n),[p1],translations=((a1,m),(a2,n)),vectors=[a1*m,a2*n],nneighbour=2)
    etime=time.time()
    print etime-stime
    a.plot(show=True)
    stime=time.time()
    b=Lattice('C'+str(m)+str(n),[p1],translations=((a1,m),(a2,n)),nneighbour=2)
    etime=time.time()
    print etime-stime
    b.plot(show=True)

def test_lattice_indices():
    p1=Point(site=0,norbital=2,nspin=2,nnambu=2,rcoord=[0.0,0.0],icoord=[0.0,0.0])
    p2=Point(site=1,norbital=2,nspin=2,nnambu=2,rcoord=[1.0,0.0],icoord=[0.0,0.0])
    a=Lattice('C',[p1,p2],nneighbour=2)
    print a.indices(nambu=True)

def test_super_lattice():
    m=20
    points=[None for i in xrange(4)]
    points[0]=Point(site=0,norbital=1,nspin=2,nnambu=1,rcoord=[0.0,0.0],icoord=[0.0,0.0])
    points[1]=Point(site=1,norbital=1,nspin=2,nnambu=1,rcoord=[0.0,1.0],icoord=[0.0,0.0])
    points[2]=Point(site=2,norbital=1,nspin=2,nnambu=1,rcoord=[1.0,0.0],icoord=[0.0,0.0])
    points[3]=Point(site=3,norbital=1,nspin=2,nnambu=1,rcoord=[1.0,1.0],icoord=[0.0,0.0])
    a1=array([2.0,0.0])
    a2=array([0.0,2.0])
    a=SuperLattice(
        name='Super',
        sublattices=[Lattice(name='sub'+'0' if i<10 else ''+str(i),points=points_shifted(points,a1*i,scope='sub'+'0' if i<10 else ''+str(i))) for i in xrange(m)],
        vectors=[a1*m,a2],
        nneighbour=2
        )
    a.plot()
    print a.indices()
    for lattice in a.sub_lattices:
        print lattice.indices()
