import sys
for arg in sys.argv:
    if arg in ('name','all'):
        from Test.Name import *
        test_name()
    if arg in ('index','all'):
        from Test.Index import *
        test_index()
    if arg in ('table','all'):
        from Test.Table import *
        test_table()
    if arg in ('basicgeometry','all'):
        from Test.BasicGeometry import *
        test_basicgeometry()
    if arg in ('indexpackage','all'):
        from Test.IndexPackage import *
        test_indexpackage()
    if arg in ('bond','all'):
        from Test.Bond import *
        test_bond()
    if arg in ('operator','all'):
        from Test.Operator import *
        test_operator()
    if arg in ('basespace','all'):
        from Test.BaseSpace import *
        test_basespace()
    if arg in ('lattice','all'):
        from Test.Lattice import *
        test_lattice()
    if arg in ('basise','all'):
        from Test.BasisE import *
        test_basise()
    if arg in ('quadratic','all'):
        from Test.Quadratic import *
        test_quadratic()
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
