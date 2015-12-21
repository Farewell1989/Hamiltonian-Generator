from Hamiltonian.Core.BasicClass.LatticePy import *
import time
def test_lattice():
    test_lattice_body() 
    test_lattice_table()
    test_super_lattice()   
    
def test_lattice_body():
    p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=Fermi(norbital=2,nspin=3,nnambu=2))
    a1=array([1.0,0.0])
    a2=array([0.0,1.0])
    m=2;n=2
    stime=time.time()
    a=Lattice('L'+str(m)+str(n),[p1],translations=((a1,m),(a2,n)),vectors=[a1*m,a2*n],nneighbour=2)
    etime=time.time()
    print etime-stime
    for p in a.points:
        print p
    for bond in a.bonds:
        print bond
    a.plot(show=True)
    stime=time.time()
    b=Lattice('C'+str(m)+str(n),[p1],translations=((a1,m),(a2,n)),nneighbour=2)
    etime=time.time()
    print etime-stime
    b.plot(show=True)

def test_lattice_table():
    p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=Fermi(norbital=2,nspin=2,nnambu=2))
    p2=Point(site=1,rcoord=[1.0,0.0],icoord=[0.0,0.0],struct=Fermi(norbital=2,nspin=2,nnambu=2))
    a=Lattice('C',[p1,p2],nneighbour=2)
    print a.table(nambu=True)

def test_super_lattice():
    m=4
    points=[None for i in xrange(4)]
    points[0]=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],struct=Fermi(norbital=1,nspin=2,nnambu=1))
    points[1]=Point(site=1,rcoord=[0.0,1.0],icoord=[0.0,0.0],struct=Fermi(norbital=1,nspin=2,nnambu=1))
    points[2]=Point(site=2,rcoord=[1.0,0.0],icoord=[0.0,0.0],struct=Fermi(norbital=1,nspin=2,nnambu=1))
    points[3]=Point(site=3,rcoord=[1.0,1.0],icoord=[0.0,0.0],struct=Fermi(norbital=1,nspin=2,nnambu=1))
    a1=array([2.0,0.0])
    a2=array([0.0,2.0])
    a=SuperLattice(
        name='Super',
        sublattices=[Lattice(name='sub'+str(i),points=translation(points,a1*i)) for i in xrange(m)],
        vectors=[a1*m,a2],
        nneighbour=2
        )
    a.plot()
    print a.table()
    for lattice in a.sublattices:
        print lattice.table()
