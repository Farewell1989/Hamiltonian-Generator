from Hamiltonian.FLQTPy import *
from BasicClass.BaseSpacePy import *
from BasicClass.LatticePy import *
def test_flqt():
    N=50
    mu1=3.0
    mu2=0.0
    name='FLQT'
    p1=Point(site=0,rcoord=[0.0],icoord=[0.0],norbital=1,nspin=1,nnambu=2,scope=name)
    a1=array([1.0])
    a=FLQT(
        name=       name,
        generator=  Generator(
            lattice=    Lattice(name=name,points=[p1],translations=[(a1,N)]),
            #lattice=    Lattice(name=name,points=[p1],vectors=[a1]),
            parameters= {'mu1':mu1,'mu2':mu2},
            terms=[     Hopping('t1',-1.0),
                        Onsite('mu',0.0,modulate=lambda **karg: mu1 if karg['t']<1 else mu2),
                        Pairing('delta',0.5,neighbour=1,amplitude=lambda bond: 1 if bond.rcoord[0]>0 else -1)
                        ],
            nambu=True
            )
        )
    #a.addapps('EB',EB(path=BaseSpace('t',array([0,1])),save_data=False,run=TBAEB))
    a.addapps('EB',EB(ts=TSpace(array([0,1,2])),save_data=False,run=FLQTEB))
    a.runapps()
