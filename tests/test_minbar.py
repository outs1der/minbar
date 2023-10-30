#!/usr/bin/env python
"""
(Py)minbar tests
"""

import os
import numpy as np
import pandas as pd
import astropy.units as u

# local modules
import minbar

import pytest

# MINBAR DR1

N_BURSTS, N_OBS, N_SRC = 7111, 118848, 115

def test_files_exist():
    """
    If this one fails, they all will
    """

    files = ['minbar.txt','minbar-obs.txt','minbar_sources.fits']

    path = os.path.join(os.path.dirname(minbar.__file__),'data')

    if not os.path.exists(path):
        raise AssertionError("MINBAR data path {} not found!".format(
            path))

    for _file in files:
        if not os.path.exists(os.path.join(path, _file)):
            raise AssertionError("""MINBAR source file {} not found!
  Check the path {}""".format(_file, path))

# This fixture is to read in the data and keep the 3 objects persistent
# for the entire test suite

@pytest.fixture(scope='session')
def load_data_fixture():

    o = minbar.Observations()
    b = o.bursts
    s = minbar.Sources()

    return (b, o, s)

# define test functions, which will be run below (if their name begins
# with "test"

def test_bursts(load_data_fixture):

    b, o, s = load_data_fixture

    # this will obviously have to be updated if we ever extend MINBAR

    b.clear()
    result = (len(b) == N_BURSTS)

    if result:
        print("PASSED")
    else:
        print("FAILED")
        raise AssertionError("MINBAR Bursts object creation issue - \
        result not as expected!")

def test_observations(load_data_fixture):

    b, o, s = load_data_fixture

    # this will obviously have to be updated if we ever extend MINBAR

    o.clear()
    result = (len(o) == N_OBS)

    if result:
        print("PASSED")
    else:
        print("FAILED")
        raise AssertionError("MINBAR Observations object creation issue - \
        result not as expected!")

def test_Sources(load_data_fixture):

    b, o, s = load_data_fixture

    s.clear()
    if len(s) != N_SRC:
        raise AssertionError("MINBAR Sources object creation issue - \
        too many/not enough sources!")

def test_selection(load_data_fixture):

    b, o, s = load_data_fixture

    b.clear()
    b.name_like('1636')
    n_1636 = len(b)

    if n_1636 != 666:
        raise AssertionError("Got unexpected number of bursts from 4U 1636-536")

    b.instr('pca')
    n_pca = len(b)

    b.PRE()
    if len(b) != 77:
        raise AssertionError("Got unexpected number of PRE bursts from 4U 1636-536")
    b.clear()

    b.name_like('1636').instr('wfc')
    n_wfc = len(b)

    b.clear()

    b.name_like('1636').instr('jemx')
    n_jemx = len(b)

    if n_1636 != n_pca+n_wfc+n_jemx:
        raise AssertionError("Bursts from 4U 1636-536 with different instruments, don't add up")

def test_online_data():

    b = minbar.Bursts()

    b.clear()
    b.name_like('1636').instr('pca').PRE()
    ids = b['entry']
    # try to get some burst data

    d = b.get_burst_data(ids[0])
    if type(d) != pd.DataFrame:
        raise AssertionError("Can't retrieve data for burst {}".format(ids[0]))

def test_distcor():

    b = minbar.Bursts()
    b.clear()

    distcor, distcor_err, dist, dist_err = b.create_distance_correction()

    b.name_like('1728')
    if (len(set(distcor[b.selection])) != 1) |\
        (not np.allclose(b.distcor[b.selection][0], 2.31641477e+45*u.cm**2)):
        raise AssertionError("Some problem with the distance corrections")


