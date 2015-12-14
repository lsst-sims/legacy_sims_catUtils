"""
SNObject_tests:
A Class containing tests to check crictical functionality for SNObject.py

The following functionality is tested:

    - SED (flambda) for unextincted SEDs in SNCosmo and SNObject
    - SED (flambda) for MW extincted SEDs in SNCosmo and SNObject (independent
        implementations of extinction using OD94 model.)
    - Band Flux for extincted SED in r Band
    - Band Mag for extincted SED in r Band

SNIaCatalog_tests:
A Class containing tests to check crictical functionality for SNIaCatalog 
"""
import os
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import unittest

# Lsst Sims Dependencies
import lsst.utils.tests as utilsTests
from lsst.sims.photUtils import Bandpass
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
import astropy
import eups

# Routines Being Tested
from lsst.sims.catUtils.mixins import SNObject
from lsst.sims.catUtils.mixins import SNIaCatalog

# External packages used
# import pandas as pd
# from pandas.util.testing import assert_frame_equal
import sncosmo
import astropy


class astropy_tests(unittest.TestCase):
    def setUp(self):
        from astropy.config import get_config_dir

        print astropy.__version__
        dir = get_config_dir() 
        print dir

    def test_config(self):
        pass


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(astropy_tests)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
