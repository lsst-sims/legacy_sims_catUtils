import unittest
import os
import numpy as np

import lsst.utils.tests as utilsTests
from lsst.utils import getPackageDir
from lsst.sims.catalogs.generation.db import fileDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.mixins import PhotometrySSM
from lsst.sims.photUtils import BandpassDict, SedList

class LSST_SSM_photCat(InstanceCatalog, PhotometrySSM):
    column_outputs = ['id', 'lsst_u' ,'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y']

    default_formats = {'f':'%.13f'}

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


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(SSMphotometryTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
