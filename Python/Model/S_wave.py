from Hamiltonian.TBAPy import *
def S_wave():
    p1=Point(site=0,rcoord=[0.0],icoord=[0.0],norbital=1,nspin=2,nnambu=2,scope='WG')
    a1=array([1.0])
    a=TBA(
        dout=       'Model/Results',
        name=       'WG',
        generator=  Generator(
            lattice=    Lattice(name='WG',points=[p1],vectors=[a1]),
            terms=[     Hopping('t1',-1.0),
                        Hopping('t2',0.0j,indexpackages=sigmay('sp'),amplitude=lambda bond: 1 if bond.rcoord[0]>0 else -1),
                        Onsite('mu',0.0,indexpackages=sigmay('sp')),
                        Pairing('delta',0.1,neighbour=0,indexpackages=sigmay('sp')*1.0j)
                        ],
            nambu=      True
            )
         )
    a.addapps('EB',EB(path=line_1d(nk=1000),save_data=False,run=TBAEB))
    a.addapps('DOS',DOS(BZ=line_1d(nk=10000),delta=0.01,ne=400,save_data=False,run=TBADOS))
    a.runapps()