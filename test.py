import sys
for arg in sys.argv:
    if arg in ('name','all'):
        from Test.Name import *
        test_name()
    if arg in ('index','all'):
        from Test.Index import *
        test_index()
        test_index_functions()
    if arg in ('table','all'):
        from Test.Table import *
        test_table()
        test_table_functions_index()
        test_table_functions_string()
    if arg in ('basicgeometry','all'):
        from Test.BasicGeometry import *
        test_point()
        test_basicgeometry_functions()
    if arg in ('indexpackage','all'):
        from Test.IndexPackage import *
        test_indexpackage()
        test_indexpackage_functions()
    if arg in ('bond','all'):
        from Test.Bond import *
        test_bond()
    if arg in ('operator','all'):
        from Test.Operator import *
        test_operator()
        test_operatorlist()
    if arg in ('basespace','all'):
        from Test.BaseSpace import *
        test_basespace()
        test_basespace_functions()
    if arg in ('lattice','all'):
        from Test.Lattice import *
        test_lattice()
    if arg in ('basise','all'):
        from Test.BasisE import *
        test_basise()
    if arg in ('quadratic','all'):
        from Test.Quadratic import *
        test_quadratic()
        test_quadratic_operators()
    if arg in ('hubbard','all'):
        from Test.Hubbard import *
        test_hubbard()
    if arg in ('optrep','all'):
        from Test.OperatorRepresentation import *
        test_opt_rep()
    if arg in ('engineapp','all'):
        from Test.EngineApp import *
        test_engineapp()
    if arg in ('lanczos','all'):
        from Test.Lanczos import *
        test_lanczos()
    if arg in ('tba','all'):
        from Test.TBA import *
        test_tba()
    if arg in ('flqt','all'):
        from Test.FLQT import *
        test_flqt()
    if arg in ('onr','all'):
        from Test.ONR import *
        test_onr()
    if arg in ('vca','all'):
        from Test.VCA import *
        test_vca()
    if arg in ('vcacct','all'):
        from Test.VCACCT import *
        test_vcacct()
    if arg in ('model.s_wave','model'):
        from Model.S_wave import *
        S_wave()
    if arg in ('model.hexagon','model'):
        from Model.Hexagon import *
        Hexagon_VCA()
        Hexagon_TBA()
