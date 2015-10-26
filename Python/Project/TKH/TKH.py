from Hamiltonian.VCAPy import *
name="TKH"
points=[None for i in xrange(12)]
points[0]=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[1]=Point(site=1,rcoord=[1.0,0.0],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[2]=Point(site=2,rcoord=[2.0,0.0],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[3]=Point(site=3,rcoord=[3.0,0.0],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[4]=Point(site=4,rcoord=[0.5,-sqrt(3.0)/2],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[5]=Point(site=5,rcoord=[1.5,-sqrt(3.0)/2],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[6]=Point(site=6,rcoord=[2.5,-sqrt(3.0)/2],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[7]=Point(site=7,rcoord=[1.0,-sqrt(3.0)],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[8]=Point(site=8,rcoord=[2.0,-sqrt(3.0)],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[9]=Point(site=9,rcoord=[0.5,sqrt(3.0)/2],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[10]=Point(site=10,rcoord=[1.5,sqrt(3.0)/2],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
points[11]=Point(site=11,rcoord=[2.5,sqrt(3.0)/2],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
a1=array([0.0,2*sqrt(3.0)])
a2=array([3.0,sqrt(3.0)])
b1=array([1.0,0.0])
b2=array([0.5,sqrt(3.0)/2])
l12=Lattice(name=name,points=points,vectors=[a1,a2])
l12_open=Lattice(name=name,points=points)
cell=Lattice(name=name,points=[points[0]],vectors=[b1,b2])

def kitaev_hopping(bond):
    theta=azimuthd(bond.rcoord)
    if abs(theta)<RZERO or abs(theta-180)<RZERO: return sigmax('sp')
    if abs(theta-60)<RZERO or abs(theta-240)<RZERO: return sigmay('sp')
    if abs(theta-120)<RZERO or abs(theta-300)<RZERO: return sigmaz('sp')

t1,t2,U=0.0,-1.0,9.0
TKH=ONR(
    name=       name,
    ensemble=   'c',
    filling=    0.5,
    mu=         U/2,
    basis=      BasisE((24,12)),
    generator=  Generator(
        lattice=    l12_open,
        terms=[     Hopping('t1',t1),
                    Hopping('t2',t2,indexpackages=kitaev_hopping),
                    Hubbard('U',[U])
                    ]
        )
    )
#TKH.addapps('GFC',GFC(nstep=100,save_data=False,run=ONRGFC))
#TKH.addapps('EB',EB(path=hexagon_gkm(nk=100),emin=-5,emax=5,ne=401,save_data=False,run=VCAEB))
#TKH.runapps()
stime=time.time()
TKH.set_matrix()
print time.time()-stime
print eigsh(TKH.matrix,k=12,which='SA',return_eigenvectors=False)
print time.time()-stime