import sys
for arg in sys.argv:
    if arg in ('name','all'):
        from BasicClass.NamePy import *
        test_name()
    if arg in ('index','all'):
        from BasicClass.IndexPy import *
        test_index()
        test_index_functions()
    if arg in ('table','all'):
        from BasicClass.TablePy import *
        test_table()
        test_table_functions_index()
        test_table_functions_string()
    if arg in ('basicgeometry','all'):
        from BasicClass.BasicGeometryPy import *
        test_point()
        test_basicgeometry_functions()
    if arg in ('indexpackage','all'):
        from BasicClass.IndexPackagePy import *
        test_indexpackage()
        test_indexpackage_functions()
    if arg in ('bond','all'):
        from BasicClass.BondPy import *
        test_bond
    if arg in ('operator','all'):
        from BasicClass.OperatorPy import *
        test_operator()
        test_operatorlist()
    if arg in ('basespace','all'):
        from BasicClass.BaseSpacePy import *
        test_basespace()
        test_basespace_functions()
    if arg in ('lattice','all'):
        from BasicClass.LatticePy import *
        test_lattice()
        test_lattice_indices()
    if arg in ('basise','all'):
        from BasicClass.BasisEPy import *
        test_basise()
    if arg in ('quadratic','all'):
        from BasicClass.QuadraticPy import *
        for i in xrange(1):
            test_quadratic()
            test_quadratic_operators()
    if arg in ('hubbard','all'):
        from BasicClass.HubbardPy import *
        test_hubbard()
    if arg in ('optrep','all'):
        from BasicClass.OperatorRepresentationPy import *
        test_opt_rep()
    if arg in ('engineapp','all'):
        from BasicClass.EngineAppPy import *
        test_engineapp()
    if arg in ('lanczos','all'):
        from BasicAlgorithm.LanczosPy import *
        test_lanczos()
    if arg in ('tba','all'):
        from Hamiltonian.TBAPy import *
        test_tba()
    if arg in ('flqt','all'):
        from Hamiltonian.FLQTPy import *
        test_flqt()
    if arg in ('onr','all'):
        from Hamiltonian.ONRPy import *
        test_onr()
    if arg in ('vca','all'):
        from Hamiltonian.VCAPy import *
        test_vca()
    if arg in ('model.s_wave','model'):
        from Model.S_wave import *
        S_wave()
    if arg in ('model.hexagon','model'):
        from Model.Hexagon import *
        #Hexagon_VCA()
        Hexagon_TBA()
    if arg in ('project.tkh'):
        from Project.TKH.Run import *