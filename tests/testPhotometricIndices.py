import os
import numpy as np
import unittest
import lsst.utils.tests as utilsTests

from lsst.utils import getPackageDir
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.generation.db import fileDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.mixins import PhotometryStars

class baselineStarCatalog(InstanceCatalog, PhotometryStars):
    column_outputs= ['raJ2000', 'decJ2000',
                    'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                    'sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r',
                    'sigma_lsst_i', 'sigma_lsst_z', 'sigma_lsst_y']

    default_formats = {'f':'%.13f'}


class uStarCatalog(InstanceCatalog, PhotometryStars):
    column_outputs = ['raJ2000', 'decJ2000', 'lsst_u', 'sigma_lsst_u']

    default_formats = {'f':'%.13f'}


class gzStarCatalog(InstanceCatalog, PhotometryStars):
    column_outputs = ['raJ2000', 'decJ2000', 'lsst_g', 'lsst_z',
                      'sigma_lsst_g', 'sigma_lsst_z']

    default_formats = {'f':'%.13f'}


class gzUncertaintyStarCatalog(InstanceCatalog, PhotometryStars):
    column_outputs = ['raJ2000', 'decJ2000', 'sigma_lsst_g', 'sigma_lsst_z']

    default_formats = {'f':'%.13f'}


class IndexTestCaseStars(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.obs = ObservationMetaData(bandpassName=['u', 'g', 'r', 'i', 'z', 'y'],
                                      m5 = [22.0, 23.0, 24.0, 25.0, 26.0, 27.0])

        baselineDtype = np.dtype([
                                 ('raJ2000', np.float),
                                 ('decJ2000', np.float),
                                 ('lsst_u', np.float),
                                 ('lsst_g', np.float),
                                 ('lsst_r', np.float),
                                 ('lsst_i', np.float),
                                 ('lsst_z', np.float),
                                 ('lsst_y', np.float),
                                 ('sigma_lsst_u', np.float),
                                 ('sigma_lsst_g', np.float),
                                 ('sigma_lsst_r', np.float),
                                 ('sigma_lsst_i', np.float),
                                 ('sigma_lsst_z', np.float),
                                 ('sigma_lsst_y', np.float),
                                 ])

        dbdtype = np.dtype([
                           ('id', np.int),
                           ('raJ2000', np.float),
                           ('decJ2000', np.float),
                           ('sedFilename', str, 100),
                           ('magNorm', np.float),
                           ('galacticAv', np.float)
                           ])

        inputDir = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'testData')
        inputFile = os.path.join(inputDir, 'IndicesTestCatalogStars.txt')

        cls.db = fileDBObject(inputFile, runtable='test',
                              idColKey='id', dtype=dbdtype)

        cat = baselineStarCatalog(cls.db, obs_metadata=cls.obs)
        cls.catName = os.path.join(getPackageDir('sims_catUtils'), 'tests',
                                   'scratchSpace', 'indicesControlCat.txt')

        cat.write_catalog(cls.catName)
        cls.controlData = np.genfromtxt(cls.catName, dtype=baselineDtype, delimiter=',')


    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.catName):
            os.unlink(cls.catName)


    def test_u_star_catalog(self):
        """
        Test that a catalog which only cares about u does not calculate any other magnitudes.
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'indicesUCat.txt')
        dtype = np.dtype([
                          ('raJ2000', np.float),
                          ('decJ2000', np.float),
                          ('lsst_u', np.float),
                          ('sigma_lsst_u', np.float)
                         ])

        cat = uStarCatalog(self.db, obs_metadata=self.obs)
        cat.write_catalog(catName)
        testData = np.genfromtxt(catName, dtype=dtype, delimiter=',')
        np.testing.assert_array_almost_equal(self.controlData['raJ2000'], testData['raJ2000'], 10)
        np.testing.assert_array_almost_equal(self.controlData['decJ2000'], testData['decJ2000'], 10)
        np.testing.assert_array_almost_equal(self.controlData['lsst_u'], testData['lsst_u'], 10)
        np.testing.assert_array_almost_equal(self.controlData['sigma_lsst_u'], testData['sigma_lsst_u'], 10)

        self.assertTrue('lsst_g' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_g' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_r' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_r' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_i' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_i' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_z' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_z' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_y' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_y' not in cat._actually_calculated_columns)

        if os.path.exists(catName):
            os.unlink(catName)


    def test_gz_star_catalog(self):
        """
        Test that a catalog which only cares about g and z does not calculate any other magnitudes
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'indicesGZCat.txt')
        dtype = np.dtype([
                          ('raJ2000', np.float),
                          ('decJ2000', np.float),
                          ('lsst_g', np.float),
                          ('lsst_z', np.float),
                          ('sigma_lsst_g', np.float),
                          ('sigma_lsst_z', np.float)
                         ])

        cat = gzStarCatalog(self.db, obs_metadata=self.obs)
        cat.write_catalog(catName)
        testData = np.genfromtxt(catName, dtype=dtype, delimiter=',')
        np.testing.assert_array_almost_equal(self.controlData['raJ2000'], testData['raJ2000'], 10)
        np.testing.assert_array_almost_equal(self.controlData['decJ2000'], testData['decJ2000'], 10)
        np.testing.assert_array_almost_equal(self.controlData['lsst_g'], testData['lsst_g'], 10)
        np.testing.assert_array_almost_equal(self.controlData['sigma_lsst_g'], testData['sigma_lsst_g'], 10)
        np.testing.assert_array_almost_equal(self.controlData['lsst_z'], testData['lsst_z'], 10)
        np.testing.assert_array_almost_equal(self.controlData['sigma_lsst_z'], testData['sigma_lsst_z'], 10)

        self.assertTrue('lsst_u' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_u' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_r' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_r' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_i' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_i' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_y' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_y' not in cat._actually_calculated_columns)

        if os.path.exists(catName):
            os.unlink(catName)


    def test_gz_uncertainty_star_catalog(self):
        """
        Test that a catalog which only cares about g and z uncertainties does not calculate any other magnitudes
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'indicesGZUncertaintyCat.txt')
        dtype = np.dtype([
                          ('raJ2000', np.float),
                          ('decJ2000', np.float),
                          ('sigma_lsst_g', np.float),
                          ('sigma_lsst_z', np.float)
                         ])

        cat = gzUncertaintyStarCatalog(self.db, obs_metadata=self.obs)
        cat.write_catalog(catName)
        testData = np.genfromtxt(catName, dtype=dtype, delimiter=',')
        np.testing.assert_array_almost_equal(self.controlData['raJ2000'], testData['raJ2000'], 10)
        np.testing.assert_array_almost_equal(self.controlData['decJ2000'], testData['decJ2000'], 10)
        np.testing.assert_array_almost_equal(self.controlData['sigma_lsst_g'], testData['sigma_lsst_g'], 10)
        np.testing.assert_array_almost_equal(self.controlData['sigma_lsst_z'], testData['sigma_lsst_z'], 10)

        self.assertTrue('lsst_u' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_u' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_r' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_r' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_i' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_i' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_y' not in cat._actually_calculated_columns)
        self.assertTrue('sigma_lsst_y' not in cat._actually_calculated_columns)
        self.assertTrue('lsst_g' in cat._actually_calculated_columns)
        self.assertTrue('lsst_z' in cat._actually_calculated_columns)

        if os.path.exists(catName):
            os.unlink(catName)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(IndexTestCaseStars)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
