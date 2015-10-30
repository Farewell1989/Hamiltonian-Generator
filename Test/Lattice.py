from BasicClass.LatticePy import *
import time
def test_lattice():
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
