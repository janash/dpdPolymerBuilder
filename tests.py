# /usr/bin/python

import makerFunctionsClass as mf
import math

sysin = '''numtypes=3
sequence=1_1_1_1-2-2-2-2
numChains=2
volFraction=0.25
angles=True
angle=1_1_1:1
angle=1_1-2:1
angle=2-2-2:2'''

x = mf.makerFunctions(sysin)

"""
Tests that parameters are correctly parsed from input file or string
"""

def test_num_types():
    assert x._numTypes == '3'

def test_num_chains():
    """
    Tests num chains assignment
    """
    assert x._numChains == '2'

def test_system_name():
    assert x._nameString == '1_1_1_1-2-2-2-2_2_0p25.data'

def test_num_types():
    assert x._nTypesinChain == 2

def test_num_bond_types():
    assert x._nBondTypes == 2

def test_bond_types():
    assert (x._nBondType == [1, 1, 1, 2, 2, 2, 2] or x._nBondType == [2, 2, 2, 1, 1, 1, 1])


"""
Tests that system parameters are calculated correctly
"""

def test_polymer_atoms():
    assert x._nM == 16

def test_box_length():
    assert x._boxLength == math.pow(16/ (3* 0.25), float(1/3))

def test_total_atoms():
    assert x._numAtoms == 64

def test_total_bonds():
    assert x._numBonds == 14
