from __future__ import with_statement
import os
import numpy as np
import unittest
import lsst.utils.tests

from lsst.utils import getPackageDir
from lsst.sims.catalogs.db import fileDBObject
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import compound
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies
from lsst.sims.photUtils import BandpassDict, Bandpass, Sed, PhotometricParameters, calcMagError_m5


def setup_module(module):
    lsst.utils.tests.init()


class PhotometryCartoon(object):

    @compound('cartoon_u', 'cartoon_g', 'cartoon_r', 'cartoon_i', 'cartoon_z')
    def get_cartoon_mags(self):
        if not hasattr(self, 'cartoonBandpassDict'):
            self.cartoonBandpassDict = BandpassDict.loadTotalBandpassesFromFiles(
                                                     bandpassNames = ['u', 'g', 'r', 'i', 'z'],
                                                     bandpassDir = os.path.join(getPackageDir('sims_photUtils'), 'tests',
                                                                               'cartoonSedTestData'),
                                                     bandpassRoot = 'test_bandpass_'
                                                     )

        return self._magnitudeGetter(self.cartoonBandpassDict, self.get_cartoon_mags._colnames)


    @compound('sigma_cartoon_u', 'sigma_cartoon_g', 'sigma_cartoon_r', 'sigma_cartoon_i', 'sigma_cartoon_z')
    def get_cartoon_uncertainty(self):
        return self._magnitudeUncertaintyGetter(['cartoon_u', 'cartoon_g', 'cartoon_r', 'cartoon_i', 'cartoon_z'],
                                                ['c_u', 'c_g', 'c_r', 'c_i', 'c_z', 'c_y'],
                                                 'cartoonBandpassDict')


class CartoonStars(InstanceCatalog, PhotometryStars, PhotometryCartoon):
    column_outputs = ['cartoon_u', 'cartoon_g',
                      'sigma_cartoon_u', 'sigma_cartoon_g']

    default_formats = {'f': '%.13f'}


class CartoonUncertaintyTestCase(unittest.TestCase):
    """
    This unit test suite will verify that our 'arbitrarily extensible' framework
    for writing magnitude uncertainty getters actually behaves as advertised
    """

    def setUp(self):
        self.normband = Bandpass()
        self.normband.imsimBandpass()
        self.uband = Bandpass()
        self.uband.readThroughput(os.path.join(getPackageDir('sims_photUtils'), 'tests', 'cartoonSedTestData', 'test_bandpass_u.dat'))
        self.gband = Bandpass()
        self.gband.readThroughput(os.path.join(getPackageDir('sims_photUtils'), 'tests', 'cartoonSedTestData', 'test_bandpass_g.dat'))

    def test_stars(self):
        obs = ObservationMetaData(bandpassName=['c_u', 'c_g'],
                                  m5=[25.0, 26.0])
        db_dtype = np.dtype([
                           ('id', np.int),
                           ('raJ2000', np.float),
                           ('decJ2000', np.float),
                           ('sedFilename', str, 100),
                           ('magNorm', np.float),
                           ('galacticAv', np.float)
                           ])

        inputDir = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'testData')
        inputFile = os.path.join(inputDir, 'IndicesTestCatalogStars.txt')
        db = fileDBObject(inputFile, dtype=db_dtype, runtable='test', idColKey='id')
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'cartoonStarCat.txt')
        cat = CartoonStars(db, obs_metadata=obs)
        cat.write_catalog(catName)

        dtype = np.dtype([(name, np.float) for name in cat.column_outputs])

        controlData = np.genfromtxt(catName, dtype=dtype, delimiter=',')

        if os.path.exists(catName):
            os.unlink(catName)

        db_columns = db.query_columns(['id', 'raJ2000', 'decJ2000', 'sedFilename', 'magNorm', 'galacticAv'])

        sedDir = os.path.join(getPackageDir('sims_sed_library'), 'starSED', 'kurucz')

        for ix, line in enumerate(db_columns.next()):
            spectrum = Sed()
            spectrum.readSED_flambda(os.path.join(sedDir, line[3]))
            fnorm = spectrum.calcFluxNorm(line[4], self.normband)
            spectrum.multiplyFluxNorm(fnorm)
            a_x, b_x = spectrum.setupCCMab()
            spectrum.addCCMDust(a_x, b_x, A_v=line[5])
            umag = spectrum.calcMag(self.uband)
            self.assertAlmostEqual(umag, controlData['cartoon_u'][ix], 3)
            gmag = spectrum.calcMag(self.gband)
            self.assertAlmostEqual(gmag, controlData['cartoon_g'][ix], 3)
            umagError, gamma = calcMagError_m5(umag, self.uband, obs.m5['c_u'], PhotometricParameters())
            gmagError, gamma = calcMagError_m5(gmag, self.gband, obs.m5['c_g'], PhotometricParameters())
            self.assertAlmostEqual(umagError, controlData['sigma_cartoon_u'][ix], 3)
            self.assertAlmostEqual(gmagError, controlData['sigma_cartoon_g'][ix], 3)


    def test_mixed_stars(self):
        """
        Here we will test the (somewhat absurd) case of a catalog with two different bandpasses
        (lsst_ and cartoon_) in order to verify that gamma values are being cached correctly
        """

        lsst_u_band = Bandpass()
        lsst_u_band.readThroughput(os.path.join(getPackageDir('throughputs'), 'baseline', 'total_u.dat'))
        lsst_g_band = Bandpass()
        lsst_g_band.readThroughput(os.path.join(getPackageDir('throughputs'), 'baseline', 'total_g.dat'))

        obs = ObservationMetaData(bandpassName=['c_u', 'c_g', 'u', 'g'],
                                  m5=[25.0, 26.0, 15.0, 16.0])
        # make the difference in m5 between the two bandpass systems extreme
        # so that, in the unit test, we can be sure that the correct values
        # are being used for the correct getters

        db_dtype = np.dtype([
                           ('id', np.int),
                           ('raJ2000', np.float),
                           ('decJ2000', np.float),
                           ('sedFilename', str, 100),
                           ('magNorm', np.float),
                           ('galacticAv', np.float)
                           ])

        inputDir = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'testData')
        inputFile = os.path.join(inputDir, 'IndicesTestCatalogStars.txt')
        db = fileDBObject(inputFile, dtype=db_dtype, runtable='test', idColKey='id')
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'cartoonStarCat.txt')
        cat = CartoonStars(db, obs_metadata=obs, column_outputs=['lsst_u', 'lsst_g', 'sigma_lsst_u', 'sigma_lsst_g'])
        cat.write_catalog(catName)

        dtype = np.dtype([(name, np.float) for name in cat._column_outputs])

        controlData = np.genfromtxt(catName, dtype=dtype, delimiter=',')

        if os.path.exists(catName):
            os.unlink(catName)

        db_columns = db.query_columns(['id', 'raJ2000', 'decJ2000', 'sedFilename', 'magNorm', 'galacticAv'])

        sedDir = os.path.join(getPackageDir('sims_sed_library'), 'starSED', 'kurucz')

        for ix, line in enumerate(db_columns.next()):
            spectrum = Sed()
            spectrum.readSED_flambda(os.path.join(sedDir, line[3]))
            fnorm = spectrum.calcFluxNorm(line[4], self.normband)
            spectrum.multiplyFluxNorm(fnorm)
            a_x, b_x = spectrum.setupCCMab()
            spectrum.addCCMDust(a_x, b_x, A_v=line[5])
            umag = spectrum.calcMag(self.uband)
            self.assertAlmostEqual(umag, controlData['cartoon_u'][ix], 3)
            gmag = spectrum.calcMag(self.gband)
            self.assertAlmostEqual(gmag, controlData['cartoon_g'][ix], 3)
            lsst_umag = spectrum.calcMag(lsst_u_band)
            self.assertAlmostEqual(lsst_umag, controlData['lsst_u'][ix], 3)
            lsst_gmag = spectrum.calcMag(lsst_g_band)
            self.assertAlmostEqual(lsst_gmag, controlData['lsst_g'][ix], 3)
            umagError, gamma = calcMagError_m5(umag, self.uband, obs.m5['c_u'], PhotometricParameters())
            gmagError, gamma = calcMagError_m5(gmag, self.gband, obs.m5['c_g'], PhotometricParameters())
            self.assertAlmostEqual(umagError, controlData['sigma_cartoon_u'][ix], 3)
            self.assertAlmostEqual(gmagError, controlData['sigma_cartoon_g'][ix], 3)

            lsst_umagError, gamma = calcMagError_m5(lsst_umag, lsst_u_band, obs.m5['u'], PhotometricParameters())
            lsst_gmagError, gamma = calcMagError_m5(lsst_gmag, lsst_g_band, obs.m5['g'], PhotometricParameters())

            self.assertAlmostEqual(lsst_umagError, controlData['sigma_lsst_u'][ix], 3)
            self.assertAlmostEqual(lsst_gmagError, controlData['sigma_lsst_g'][ix], 3)
            self.assertGreater(np.abs(lsst_umagError-umagError), 0.01)
            self.assertGreater(np.abs(lsst_gmagError-gmagError), 0.01)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
