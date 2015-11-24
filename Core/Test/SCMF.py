from Hamiltonian.Core.CoreAlgorithm.SCMFPy import *
from Hamiltonian.Core.BasicClass.LatticePy import *
from Hamiltonian.Core.BasicClass.BaseSpacePy import *
from Hamiltonian.DataBase.Hexagon import *
def test_scmf():
    H2=HexagonDataBase('H2',norbital=1,nspin=2,nnambu=1)
    h2=SCMF(
        name=       'H2_SCMF',
        lattice=    Lattice(name='H2',points=H2.points,vectors=H2.vectors,nneighbour=1),
        mu=         0,
        filling=    0.5,
        terms=      [
                    Hopping('t1',-1.0),
                    ],
        orders=     [
                    Onsite('afm',1.0,indexpackages=sigmaz('sp')*sigmaz('sl'),modulate=lambda **karg: karg['afm'] if 'afm' in karg else None)
                    ],
        nambu=      False
        )
    print h2.op_m
    #h2.addapps('EB',EB(hexagon_gkm(nk=100),save_data=False,plot=True,show=True,run=TBAEB))
    #h2.addapps('CP',CP(KSpace(reciprocals=h2.lattice.reciprocals,nk=100),run=TBACP))
    #h2.addapps('CN',CN(KSpace(reciprocals=h2.lattice.reciprocals,nk=200),d=10**-6,save_data=False,plot=True,show=True,run=TBACN))
    #h2.runapps()
