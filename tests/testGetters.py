import os
import numpy
import unittest
import eups

import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB
from lsst.sims.catalogs.measures.instance import InstanceCatalog, defaultSpecMap
from lsst.sims.catUtils.utils import testStarsDBObj
from lsst.sims.photUtils import Sed, Bandpass, PhotometricDefaults
from lsst.sims.photUtils import PhotometryStars

class testStarCatalog(InstanceCatalog, PhotometryStars):
    sig2sys = 0.0003

    column_outputs = ['raJ2000', 'decJ2000',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r', 'sigma_lsst_i',
                      'sigma_lsst_z', 'sigma_lsst_y', 'sedFilename', 'magNorm']

class testPhotometricUncertaintyGetters(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dbName = 'uncertaintyTestDB.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        cls.obs_metadata = makePhoSimTestDB(filename=cls.dbName, size=100, radius = 5.0)
        m5 = {'u':23.0, 'g':24.0, 'r':21.0, 'i':22.3, 'z':23.7, 'y':24.5}
        cls.obs_metadata.m5value = m5
        cls.connectionString = 'sqlite:///'+cls.dbName

        cls.skySeds = []
        cls.hardwareBandpasses = []
        cls.totalBandpasses = []
        cls.bandpasses = ['u', 'g', 'r', 'i', 'z', 'y']

        components = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                      'lens1.dat', 'lens2.dat', 'lens3.dat']

        for b in cls.bandpasses:
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughput(os.path.join(eups.productDir('throughputs'),
                                                      'baseline', 'total_%s.dat' % b))
            cls.totalBandpasses.append(bandpassDummy)

        for b in cls.bandpasses:
            finalComponents = []
            for c in components:
                finalComponents.append(os.path.join(eups.productDir('throughputs'), 'baseline', c))
            finalComponents.append(os.path.join(eups.productDir('throughputs'), 'baseline', 'filter_%s.dat' %b))
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(finalComponents)
            cls.hardwareBandpasses.append(bandpassDummy)

        for i in range(len(cls.bandpasses)):
            sedDummy = Sed()
            sedDummy.readSED_flambda(os.path.join(eups.productDir('throughputs'), 'baseline', 'darksky.dat'))
            cls.totalBandpasses[i].setM5(cls.obs_metadata.m5(cls.bandpasses[i]), sedDummy,
                                         cls.hardwareBandpasses[i],
                                         seeing=PhotometricDefaults.seeing[cls.bandpasses[i]])
            cls.skySeds.append(sedDummy)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)
        del cls.dbName
        del cls.connectionString
        del cls.obs_metadata

    def testStellarPhotometricUncertainties(self):
        starDB = testStarsDBObj(address=self.connectionString)
        starCat = testStarCatalog(starDB, obs_metadata=self.obs_metadata)

        ct = 0
        for line in starCat.iter_catalog():
            starSed = Sed()
            starSed.readSED_flambda(os.path.join(eups.productDir('sims_sed_library'),
                                                 defaultSpecMap[line[14]]))
            imsimband = Bandpass()
            imsimband.imsimBandpass()
            fNorm = starSed.calcFluxNorm(line[15], imsimband)
            starSed.multiplyFluxNorm(fNorm)

            for i in range(len(self.bandpasses)):
                controlSNR = starSed.calcSNR_psf(self.totalBandpasses[i],
                                                 self.skySeds[i],
                                                 self.hardwareBandpasses[i],
                                                 seeing=PhotometricDefaults.seeing[self.bandpasses[i]])
                controlNoverS = 1.0/(controlSNR*controlSNR) + starCat.sig2sys
                controlSigma = 2.5*numpy.log10(1.0+numpy.sqrt(controlNoverS))
                testSigma = line[8+i]
                self.assertAlmostEqual(controlSigma, testSigma, 10)
                ct += 1
        self.assertTrue(ct>0)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(testPhotometricUncertaintyGetters)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
