from Hamiltonian.TBAPy import *
def S_wave():
    p1=Point(site=0,rcoord=[0.0],icoord=[0.0],norbital=1,nspin=2,nnambu=2,scope='WG')
    a1=array([1.0])
    m=10
    a=TBA(
        dout=       'Results',
        name=       'WG',
        lattice=    Lattice(name='WG',points=[p1],vectors=[a1]),
#        lattice=    Lattice(name='WG',points=[p1],translations=[(a1,m)]),
        hopping=    [Hopping(-1.0,neighbour=1),Hopping(0.0j,neighbour=1,indexpackages=sigmay('sp'),amplitude=lambda bond: 1 if bond.rcoord[0]>0 else -1)],
        onsite=     [Onsite(0.0,indexpackages=sigmay('sp'))],
        pairing=    [Pairing(0.0,neighbour=0,indexpackages=sigmay('sp')*1.0j)]
        )
    a.get_ready()
#    print a.operators['hpst']
#    print a.matrix()
#    a.addapps('EB',EB(run=TBAEB))
    a.addapps('EB',EB(path=line_1d(nk=1000),run=TBAEB))
    a.addapps('DOS',DOS(BZ=line_1d(nk=10000),delta=0.01,ne=400,run=TBADOS))
    a.runapps()