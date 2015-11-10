import sys
for arg in sys.argv:
    if arg in ('name','all'):
        from Hamiltonian.Core.Test.Name import *
        test_name()
    if arg in ('index','all'):
        from Hamiltonian.Core.Test.Index import *
        test_index()
    if arg in ('table','all'):
        from Hamiltonian.Core.Test.Table import *
        test_table()
    if arg in ('basicgeometry','all'):
        from Hamiltonian.Core.Test.BasicGeometry import *
        test_basicgeometry()
    if arg in ('indexpackage','all'):
        from Hamiltonian.Core.Test.IndexPackage import *
        test_indexpackage()
    if arg in ('bond','all'):
        from Hamiltonian.Core.Test.Bond import *
        test_bond()
    if arg in ('operator','all'):
        from Hamiltonian.Core.Test.Operator import *
        test_operator()
    if arg in ('basespace','all'):
        from Hamiltonian.Core.Test.BaseSpace import *
        test_basespace()
    if arg in ('lattice','all'):
        from Hamiltonian.Core.Test.Lattice import *
        test_lattice()
    if arg in ('basise','all'):
        from Hamiltonian.Core.Test.BasisE import *
        test_basise()
    if arg in ('quadratic','all'):
        from Hamiltonian.Core.Test.Quadratic import *
        test_quadratic()
    if arg in ('hubbard','all'):
        from Hamiltonian.Core.Test.Hubbard import *
        test_hubbard()
    if arg in ('optrep','all'):
        from Hamiltonian.Core.Test.OperatorRepresentation import *
        test_opt_rep()
    if arg in ('engineapp','all'):
        from Hamiltonian.Core.Test.EngineApp import *
        test_engineapp()
    if arg in ('lanczos','all'):
        from Hamiltonian.Core.Test.Lanczos import *
        test_lanczos()
    if arg in ('tba','all'):
        from Hamiltonian.Core.Test.TBA import *
        test_tba()
    if arg in ('flqt','all'):
        from Hamiltonian.Core.Test.FLQT import *
        test_flqt()
    if arg in ('onr','all'):
        from Hamiltonian.Core.Test.ONR import *
        test_onr()
    if arg in ('vca','all'):
        from Hamiltonian.Core.Test.VCA import *
        test_vca()
    if arg in ('vcacct','all'):
        from Hamiltonian.Core.Test.VCACCT import *
        test_vcacct()
