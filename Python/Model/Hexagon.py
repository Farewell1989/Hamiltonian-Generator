from Hamiltonian.VCAPy import *
def Hexagon():
    t,U=-1.0,2.0
    name='Hexagon'
    m,n=1,1
    p1=Point(site=0,rcoord=[0.0,0.0],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
    p2=Point(site=1,rcoord=[0.0,sqrt(3)/3],icoord=[0.0,0.0],norbital=1,nspin=2,nnambu=1,scope=name)
    a1=array([1.0,0.0])
    a2=array([0.5,sqrt(3)/2])
    a=VCA(
        din=        'Model/Results/Coeff',
        dout=       'Model/Results',
        name=       name,
        ensemble=   'c',
        filling=    0.5,
        mu=         U/2,
        basis=      BasisE(up=(2*m*n,m*n),down=(2*m*n,m*n)),
        nspin=      1,
        cell=       Lattice(name=name,points=[p1,p2],vectors=[a1,a2]),
        generator=  Generator(
            lattice=    Lattice(name=name,points=[p1,p2],translations=[(a1,m),(a2,n)],vectors=[a1*m,a2*n]),
            terms=[     Hopping('t',t),
                        Hubbard('U',[U])
                        ]
            )
        )
    a.addapps('GFC',GFC(nstep=200,save_data=False,vtype='SY',run=ONRGFC))
    a.addapps('EB',EB(path=hexagon_gkm(nk=100),emin=-4,emax=4,ne=400,run=VCAEB,save_data=False))
    a.runapps()
    