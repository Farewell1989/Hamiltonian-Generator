from Hamiltonian.ONRPy import *
from BasicClass.LatticePy import *
def test_onr():
    U=2.0
    t=-1.0
    m=2;n=4
    p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],atom=0,norbital=1,nspin=2,nnambu=1,scope='WG'+str(m)+str(n))
    a1=array([1.0,0.0])
    a2=array([0.0,1.0])
    a=ONR(
            name=       'WG'+str(m)+str(n),
            ensemble=   'c',
            filling=    0.5,
            mu=         U/2,
            basis=      BasisE(up=(m*n,m*n/2),down=(m*n,m*n/2)),
            #basis=      BasisE((2*m*n,m*n)),
            nspin=      2,
            generator=  Generator(
                lattice=    Lattice(name='WG'+str(m)+str(n),points=[p1],translations=[(a1,m),(a2,n)]),
                terms=[     Hopping('t',t,neighbour=1),
                            Hubbard('U',[U])
                            ]
            )
        )
    a.addapps('GFC',GFC(nstep=200,save_data=False,vtype='SY',run=ONRGFC))
    a.addapps('DOS',DOS(emin=-5,emax=5,ne=401,delta=0.05,save_data=False,run=ONRDOS,show=True))
    a.runapps()
