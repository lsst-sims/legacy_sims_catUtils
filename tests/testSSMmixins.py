from __future__ import with_statement
import unittest
import os
import numpy as np
import lsst

import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.db import fileDBObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import compound
from lsst.sims.catUtils.mixins import PhotometrySSM, AstrometrySSM
from lsst.sims.photUtils import BandpassDict, SedList, PhotometricParameters
from lsst.sims.utils import _observedFromICRS

class LSST_SSM_photCat(InstanceCatalog, PhotometrySSM):
    column_outputs = ['id', 'lsst_u' ,'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y']

    default_formats = {'f':'%.13f'}


class Compound_SSM_photCat(InstanceCatalog, PhotometrySSM):
    column_outputs = ['id', 'lsst_u' ,'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'cartoon_u', 'cartoon_g', 'cartoon_r', 'cartoon_i', 'cartoon_z']

    default_formats = {'f':'%.13f'}

    @compound('cartoon_u', 'cartoon_g', 'cartoon_r', 'cartoon_i', 'cartoon_z')
    def get_cartoon_mags(self):

        if not hasattr(self, 'cartoonBandpassDict'):
            bandpassDir = os.path.join(getPackageDir('sims_photUtils'), 'tests', 'cartoonSedTestData')

            self.cartoonBandpassDict = BandpassDict.loadTotalBandpassesFromFiles(
                                                      ['u', 'g', 'r', 'i', 'z'],
                                                      bandpassDir = bandpassDir,
                                                      bandpassRoot = 'test_bandpass_'
                                                      )


        return self._magnitudeGetter(self.cartoonBandpassDict, self.get_cartoon_mags._colnames,
                                     bandpassTag='cartoon')


class SSM_dmagCat(InstanceCatalog, PhotometrySSM):
    column_outputs = ['id', 'dmagTrailing', 'dmagDetection']

    default_formats = {'f': '%.13f'}


class SSMphotometryTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dbFile = os.path.join(getPackageDir('sims_catUtils'),
                        'tests', 'testData', 'SSMphotometryCatalog.txt')

        cls.dtype = np.dtype([('id', np.int), ('sedFilename', str, 100), ('magNorm', np.float),
                              ('velRa', np.float), ('velDec', np.float)])

        cls.photDB = fileDBObject(cls.dbFile, runtable='test', dtype=cls.dtype, idColKey='id')


    def testLSSTmags(self):
        """
        Test that PhotometrySSM properly calculates LSST magnitudes
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'lsstSsmPhotCat.txt')

        cat=LSST_SSM_photCat(self.photDB)
        cat.write_catalog(catName)

        dtype = np.dtype([('id', np.int), ('u', np.float), ('g', np.float),
                          ('r', np.float), ('i', np.float), ('z', np.float),
                          ('y', np.float)])

        testData = np.genfromtxt(catName, dtype=dtype, delimiter=',')
        self.assertGreater(len(testData), 0)

        controlData = np.genfromtxt(self.dbFile, dtype=self.dtype)
        self.assertGreater(len(controlData), 0)

        LSSTbandpasses = BandpassDict.loadTotalBandpassesFromFiles()
        controlSedList = SedList(controlData['sedFilename'], controlData['magNorm'],
                                 wavelenMatch=LSSTbandpasses.wavelenMatch)

        controlMags = LSSTbandpasses.magListForSedList(controlSedList)

        for ii in range(len(controlMags)):
            for jj, bpName in enumerate(['u', 'g', 'r', 'i', 'z', 'y']):
                self.assertAlmostEqual(controlMags[ii][jj], testData[bpName][ii], 10)

        if os.path.exists(catName):
            os.unlink(catName)


    def testManyMagSystems(self):
        """
        Test that the SSM photometry mixin can simultaneously calculate magnitudes
        in multiple bandpass systems
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'compoundSsmPhotCat.txt')

        cat=Compound_SSM_photCat(self.photDB)
        cat.write_catalog(catName)

        dtype = np.dtype([('id', np.int), ('lsst_u', np.float), ('lsst_g', np.float),
                          ('lsst_r', np.float), ('lsst_i', np.float), ('lsst_z', np.float),
                          ('lsst_y', np.float),
                          ('cartoon_u', np.float), ('cartoon_g', np.float),
                          ('cartoon_r', np.float), ('cartoon_i', np.float),
                          ('cartoon_z', np.float)])

        testData = np.genfromtxt(catName, dtype=dtype, delimiter=',')
        self.assertGreater(len(testData), 0)

        controlData = np.genfromtxt(self.dbFile, dtype=self.dtype)
        self.assertGreater(len(controlData), 0)

        LSSTbandpasses = BandpassDict.loadTotalBandpassesFromFiles()
        cartoonBandpasses = BandpassDict.loadTotalBandpassesFromFiles(
                                  ['u', 'g', 'r', 'i', 'z'],
                                  bandpassDir = os.path.join(getPackageDir('sims_photUtils'), 'tests', 'cartoonSedTestData'),
                                  bandpassRoot = 'test_bandpass_'
                                  )

        controlSedList = SedList(controlData['sedFilename'], controlData['magNorm'],
                                 wavelenMatch=LSSTbandpasses.wavelenMatch)

        controlLsstMags = LSSTbandpasses.magListForSedList(controlSedList)
        controlCartoonMags = cartoonBandpasses.magListForSedList(controlSedList)

        for ii in range(len(controlLsstMags)):
            for jj, bpName in enumerate(['lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y']):
                self.assertAlmostEqual(controlLsstMags[ii][jj], testData[bpName][ii], 10)
            for jj, bpName in enumerate(['cartoon_u', 'cartoon_g', 'cartoon_r', 'cartoon_i', 'cartoon_z']):
                self.assertAlmostEqual(controlCartoonMags[ii][jj], testData[bpName][ii], 10)

        if os.path.exists(catName):
            os.unlink(catName)


    def testDmagExceptions(self):
        """
        Test that the dmagTrailing and dmagDetection getters raise expected
        exceptions
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'ssmDmagCatExceptions.txt')


        obs = ObservationMetaData()
        with self.assertRaises(RuntimeError) as context:
            cat = SSM_dmagCat(self.photDB, obs_metadata=obs)
        self.assertIn("does not specify seeing", context.exception.args[0])

        obs = ObservationMetaData(bandpassName = ['u', 'g'], seeing=[0.6, 0.5])

        with self.assertRaises(RuntimeError) as context:
            cat = SSM_dmagCat(self.photDB, obs_metadata=obs)
        self.assertIn("multiple seeing values", context.exception.args[0])

        obs = ObservationMetaData(bandpassName = 'u', seeing=0.7)
        with self.assertRaises(RuntimeError) as context:
            cat = SSM_dmagCat(self.photDB, obs_metadata=obs)
            cat.photParams = None
            cat.write_catalog(catName)
        self.assertIn("does not have an associated PhotometricParameters",
                       context.exception.args[0])

        if os.path.exists(catName):
            os.unlink(catName)


    def testDmag(self):
        """
        Test the calculation of dmagTrailing and dmagDetection
        """

        obs = ObservationMetaData(bandpassName = 'u', seeing=1.48)
        photParams = PhotometricParameters()
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'ssmDmagCat.txt')

        controlData = np.genfromtxt(self.dbFile, dtype=self.dtype)

        cat = SSM_dmagCat(self.photDB, obs_metadata=obs)
        cat.write_catalog(catName)

        dtype = np.dtype([('id', np.int), ('dmagTrail', np.float), ('dmagDetect', np.float)])
        testData = np.genfromtxt(catName, dtype=dtype, delimiter=',')
        self.assertGreater(len(testData), 0)

        a_trail = 0.76
        b_trail = 1.16
        a_det = 0.42
        b_det = 0.00

        velocity = np.sqrt(np.power(np.degrees(controlData['velRa']), 2) + np.power(np.degrees(controlData['velDec']), 2))
        x = velocity * photParams.nexp * photParams.exptime /(obs.seeing[obs.bandpass] * 24.0)
        xsq = np.power(x, 2)

        dmagTrailControl = 1.25*np.log10(1.0 + a_trail*xsq/(1.0+b_trail*x))
        dmagDetectControl = 1.25*np.log10(1.0 + a_det*xsq/(1.0+b_det*x))

        # Check against precalculated numbers, just to verify units/conversions, etc.
        self.assertTrue(abs(dmagTrailControl[-1] - 0.0219) < 0.01)
        self.assertTrue(abs(dmagDetectControl[-1] - 0.01594) < 0.01)

        np.testing.assert_array_almost_equal(dmagTrailControl, testData['dmagTrail'], 10)
        np.testing.assert_array_almost_equal(dmagDetectControl, testData['dmagDetect'], 10)

        if os.path.exists(catName):
            os.unlink(catName)


class SSM_astrometryCat(InstanceCatalog, AstrometrySSM):
    column_outputs = ['id', 'raObserved', 'decObserved']

    default_formats = {'f': '%.13f'}


class SSM_velocityCat(InstanceCatalog, AstrometrySSM):
    column_outputs = ['id', 'skyVelocity']

    default_formats= {'f': '%.13f'}


class SSMastrometryTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dbFile = os.path.join(getPackageDir('sims_catUtils'),
                        'tests', 'testData', 'SSMastrometryCatalog.txt')

        cls.dtype = np.dtype([('id', np.int), ('raJ2000', np.float), ('decJ2000', np.float),
                              ('velRa', np.float), ('velDec', np.float)])

        cls.astDB = fileDBObject(cls.dbFile, runtable='test', dtype=cls.dtype, idColKey='id')


    def testObservedRaDec(self):
        """
        Test that the mixins provided in Astrometry SSM really do convert ICRS RA, Dec
        into observed RA, Dec
        """

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'ssmAstrometryCat.txt')

        dtype = np.dtype([('id', np.int),
                          ('raObserved', np.float), ('decObserved', np.float)])

        controlData = np.genfromtxt(self.dbFile, dtype=self.dtype)

        np.random.seed(42)
        nTests = 5
        raList = np.random.random_sample(nTests)*2.0*np.pi
        decList = (np.random.random_sample(nTests)-0.5)*np.pi
        mjdList = np.random.random_sample(nTests)*5000.0 + 53850.0
        for raPointing, decPointing, mjd in zip(raList, decList, mjdList):
            obs = ObservationMetaData(pointingRA=raPointing, pointingDec=decPointing, mjd=mjd)

            cat = SSM_astrometryCat(self.astDB, obs_metadata=obs)
            cat.write_catalog(catName)

            testData = np.genfromtxt(catName, dtype=dtype, delimiter=',')
            self.assertGreater(len(testData), 0)

            raObservedControl, decObservedControl = _observedFromICRS(controlData['raJ2000'], controlData['decJ2000'],
                                                                      obs_metadata=obs, epoch=2000.0,
                                                                      includeRefraction=True)

            np.testing.assert_array_almost_equal(raObservedControl, testData['raObserved'], 10)
            np.testing.assert_array_almost_equal(decObservedControl, testData['decObserved'], 10)

            if os.path.exists(catName):
                os.unlink(catName)


    def testSkyVelocity(self):
        """
        Test that getter for sky velocity correctly calculates its output
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'ssmVelocityCat.txt')

        controlData = np.genfromtxt(self.dbFile, dtype=self.dtype)
        controlVel = np.sqrt(np.power(controlData['velRa'], 2) + np.power(controlData['velDec'], 2))

        cat = SSM_velocityCat(self.astDB)
        cat.write_catalog(catName)
        dtype = np.dtype([('id', np.int), ('vel', np.float)])
        testData = np.genfromtxt(catName, dtype=dtype)
        self.assertGreater(len(testData), 0)

        np.testing.assert_array_almost_equal(testData['vel'], controlVel, 10)

        if os.path.exists(catName):
            os.unlink(catName)

class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
