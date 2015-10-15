import unittest
import os
import numpy as np

import lsst.utils.tests as utilsTests
from lsst.utils import getPackageDir
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.generation.db import fileDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.catUtils.mixins import PhotometrySSM, AstrometrySSM
from lsst.sims.photUtils import BandpassDict, SedList
from lsst.sims.coordUtils import _observedFromICRS

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


        return self._magnitudeGetter(self.cartoonBandpassDict, bandpassTag='cartoon')


class SSMphotometryTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dbFile = os.path.join(getPackageDir('sims_catUtils'),
                        'tests', 'testData', 'SSMphotometryCatalog.txt')

        cls.dtype = np.dtype([('id', np.int), ('sedFilename', str, 100), ('magNorm', np.float)])

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

        controlData = np.genfromtxt(self.dbFile, dtype=self.dtype)

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

        controlData = np.genfromtxt(self.dbFile, dtype=self.dtype)

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



class SSM_astrometryCat(InstanceCatalog, AstrometrySSM):
    column_outputs = ['id', 'raObserved', 'decObserved', 'raPhoSim', 'decPhoSim']

    default_formats = {'f': '%.13f'}

class SSMastrometryTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dbFile = os.path.join(getPackageDir('sims_catUtils'),
                        'tests', 'testData', 'SSMastrometryCatalog.txt')

        cls.dtype = np.dtype([('id', np.int), ('raJ2000', np.float), ('decJ2000', np.float)])

        cls.astDB = fileDBObject(cls.dbFile, runtable='test', dtype=cls.dtype, idColKey='id')


    def testObservedRaDec(self):
        """
        Test that the mixins provided in Astrometry SSM really do convert ICRS RA, Dec
        into observed RA, Dec
        """

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'ssmAstrometryCat.txt')

        dtype = np.dtype([('id', np.int),
                          ('raObserved', np.float), ('decObserved', np.float),
                          ('raPhoSim', np.float), ('decPhoSim', np.float)])

        controlData = np.genfromtxt(self.dbFile, dtype=self.dtype)

        np.random.seed(42)
        nTests = 5
        raList = np.random.random_sample(nTests)*2.0*np.pi
        decList = (np.random.random_sample(nTests)-0.5)*np.pi
        mjdList = np.random.random_sample(nTests)*5000.0 + 53850.0
        for raPointing, decPointing, mjd in zip(raList, decList, mjdList):
            obs = ObservationMetaData(unrefractedRA=raPointing, unrefractedDec=decPointing, mjd=mjd)

            cat = SSM_astrometryCat(self.astDB, obs_metadata=obs)
            cat.write_catalog(catName)

            testData = np.genfromtxt(catName, dtype=dtype, delimiter=',')
            raPhoSimControl, decPhoSimControl = _observedFromICRS(controlData['raJ2000'], controlData['decJ2000'],
                                                                  obs_metadata=obs, epoch=2000.0,
                                                                  includeRefraction=False)

            np.testing.assert_array_almost_equal(raPhoSimControl, testData['raPhoSim'], 10)
            np.testing.assert_array_almost_equal(decPhoSimControl, testData['decPhoSim'], 10)

            raObservedControl, decObservedControl = _observedFromICRS(controlData['raJ2000'], controlData['decJ2000'],
                                                                      obs_metadata=obs, epoch=2000.0,
                                                                      includeRefraction=True)

            np.testing.assert_array_almost_equal(raObservedControl, testData['raObserved'], 10)
            np.testing.assert_array_almost_equal(decObservedControl, testData['decObserved'], 10)

            if os.path.exists(catName):
                os.unlink(catName)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(SSMphotometryTest)
    suites += unittest.makeSuite(SSMastrometryTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
