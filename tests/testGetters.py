import os
import numpy
import unittest

import lsst.utils
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.utils import defaultSpecMap
from lsst.sims.catUtils.utils import testStarsDBObj, testGalaxyTileDBObj
from lsst.sims.photUtils import Sed, Bandpass, LSSTdefaults, calcMagError_sed
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils.utils import setM5
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies

class testStarCatalog(InstanceCatalog, PhotometryStars):

    column_outputs = ['raJ2000', 'decJ2000',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r', 'sigma_lsst_i',
                      'sigma_lsst_z', 'sigma_lsst_y', 'sedFilename', 'magNorm', 'galacticAv']

class testGalaxyCatalog(InstanceCatalog, PhotometryGalaxies):

    column_outputs = ['raJ2000', 'decJ2000',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
                      'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
                      'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn',
                      'sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r', 'sigma_lsst_i',
                      'sigma_lsst_z', 'sigma_lsst_y',
                      'sigma_uBulge', 'sigma_gBulge', 'sigma_rBulge', 'sigma_iBulge',
                      'sigma_zBulge', 'sigma_yBulge',
                      'sigma_uDisk', 'sigma_gDisk', 'sigma_rDisk', 'sigma_iDisk',
                      'sigma_zDisk', 'sigma_yDisk',
                      'sigma_uAgn', 'sigma_gAgn', 'sigma_rAgn', 'sigma_iAgn',
                      'sigma_zAgn', 'sigma_yAgn',
                      'sedFilenameBulge', 'sedFilenameDisk', 'sedFilenameAgn',
                      'magNormBulge', 'magNormDisk', 'magNormAgn',
                      'internalAvBulge', 'internalAvDisk', 'redshift']

class testPhotometricUncertaintyGetters(unittest.TestCase):
    """
    This class will test that the getters for photometric uncertainties
    are calculating the correct values.

    This test is here rather than in sims_photUtils because the infrastructure
    to generate fake databases is in sims_catUtils.
    """

    @classmethod
    def setUpClass(cls):
        lsstDefaults=LSSTdefaults()
        cls.dbName = 'uncertaintyTestDB.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        default_obs_metadata = makePhoSimTestDB(filename=cls.dbName, size=10, radius = 5.0)
        bandpass = ['u', 'g', 'r', 'i', 'z', 'y']
        m5 = lsstDefaults._m5.values()

        cls.obs_metadata = ObservationMetaData(
                                              pointingRA = default_obs_metadata.pointingRA,
                                              pointingDec = default_obs_metadata.pointingDec,
                                              rotSkyPos = default_obs_metadata.rotSkyPos,
                                              bandpassName = bandpass,
                                              m5 = m5
                                              )

        cls.obs_metadata.setBandpassM5andSeeing(bandpassName=bandpass, m5=m5)
        cls.driver = 'sqlite'
        cls.host = ''

        cls.skySeds = []
        cls.hardwareBandpasses = []
        cls.totalBandpasses = []
        cls.bandpasses = ['u', 'g', 'r', 'i', 'z', 'y']

        components = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                      'lens1.dat', 'lens2.dat', 'lens3.dat']

        for b in cls.bandpasses:
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughput(os.path.join(lsst.utils.getPackageDir('throughputs'),
                                                      'baseline', 'total_%s.dat' % b))
            cls.totalBandpasses.append(bandpassDummy)

        for b in cls.bandpasses:
            finalComponents = []
            for c in components:
                finalComponents.append(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', c))
            finalComponents.append(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'filter_%s.dat' %b))
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(finalComponents)
            cls.hardwareBandpasses.append(bandpassDummy)

        for i in range(len(cls.bandpasses)):
            sedDummy = Sed()
            sedDummy.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSedDummy = setM5(cls.obs_metadata.m5[cls.bandpasses[i]], sedDummy,
                                       cls.totalBandpasses[i], cls.hardwareBandpasses[i],
                                       seeing=lsstDefaults.seeing(cls.bandpasses[i]),
                                       photParams=PhotometricParameters())

            cls.skySeds.append(normalizedSedDummy)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)
        del cls.dbName
        del cls.driver
        del cls.host
        del cls.obs_metadata
        del cls.totalBandpasses
        del cls.hardwareBandpasses
        del cls.skySeds

    def testStellarPhotometricUncertainties(self):
        """
        Test in the case of a catalog of stars
        """
        lsstDefaults = LSSTdefaults()
        starDB = testStarsDBObj(driver=self.driver, host=self.host, database=self.dbName)
        starCat = testStarCatalog(starDB, obs_metadata=self.obs_metadata)
        phot = PhotometryStars()

        ct = 0
        for line in starCat.iter_catalog():
            starSed = Sed()
            starSed.readSED_flambda(os.path.join(lsst.utils.getPackageDir('sims_sed_library'),
                                                 defaultSpecMap[line[14]]))
            imsimband = Bandpass()
            imsimband.imsimBandpass()
            fNorm = starSed.calcFluxNorm(line[15], imsimband)
            starSed.multiplyFluxNorm(fNorm)

            aV = numpy.float(line[16])
            a_int, b_int = starSed.setupCCMab()
            starSed.addCCMDust(a_int, b_int, A_v=aV)

            for i in range(len(self.bandpasses)):
                controlSigma = calcMagError_sed(starSed, self.totalBandpasses[i],
                                             self.skySeds[i],
                                             self.hardwareBandpasses[i],
                                             seeing=lsstDefaults.seeing(self.bandpasses[i]),
                                             photParams=PhotometricParameters())

                testSigma = line[8+i]
                self.assertAlmostEqual(controlSigma, testSigma, 10)
                ct += 1
        self.assertGreater(ct, 0)

    def testGalaxyPhotometricUncertainties(self):
        """
        Test in the case of a catalog of galaxies
        """
        lsstDefaults = LSSTdefaults()
        phot = PhotometryGalaxies()
        galDB = testGalaxyTileDBObj(driver=self.driver, host=self.host, database=self.dbName)
        galCat = testGalaxyCatalog(galDB, obs_metadata=self.obs_metadata)
        imsimband = Bandpass()
        imsimband.imsimBandpass()
        ct = 0
        for line in galCat.iter_catalog():
            bulgeSedName = line[50]
            diskSedName = line[51]
            agnSedName = line[52]
            magNormBulge = line[53]
            magNormDisk = line[54]
            magNormAgn = line[55]
            avBulge = line[56]
            avDisk = line[57]
            redshift = line[58]

            bulgeSed = Sed()
            bulgeSed.readSED_flambda(os.path.join(lsst.utils.getPackageDir('sims_sed_library'),
                                     defaultSpecMap[bulgeSedName]))
            fNorm=bulgeSed.calcFluxNorm(magNormBulge, imsimband)
            bulgeSed.multiplyFluxNorm(fNorm)

            diskSed = Sed()
            diskSed.readSED_flambda(os.path.join(lsst.utils.getPackageDir('sims_sed_library'),
                                    defaultSpecMap[diskSedName]))
            fNorm = diskSed.calcFluxNorm(magNormDisk, imsimband)
            diskSed.multiplyFluxNorm(fNorm)

            agnSed = Sed()
            agnSed.readSED_flambda(os.path.join(lsst.utils.getPackageDir('sims_sed_library'),
                                   defaultSpecMap[agnSedName]))
            fNorm = agnSed.calcFluxNorm(magNormAgn, imsimband)
            agnSed.multiplyFluxNorm(fNorm)

            a_int, b_int = bulgeSed.setupCCMab()
            bulgeSed.addCCMDust(a_int, b_int, A_v=avBulge)

            a_int, b_int = diskSed.setupCCMab()
            diskSed.addCCMDust(a_int, b_int, A_v=avDisk)

            bulgeSed.redshiftSED(redshift, dimming=True)
            diskSed.redshiftSED(redshift, dimming=True)
            agnSed.redshiftSED(redshift, dimming=True)

            bulgeSed.resampleSED(wavelen_match=self.totalBandpasses[0].wavelen)
            diskSed.resampleSED(wavelen_match=bulgeSed.wavelen)
            agnSed.resampleSED(wavelen_match=bulgeSed.wavelen)

            numpy.testing.assert_almost_equal(bulgeSed.wavelen, diskSed.wavelen)
            numpy.testing.assert_almost_equal(bulgeSed.wavelen, agnSed.wavelen)

            fl = bulgeSed.flambda + diskSed.flambda + agnSed.flambda

            totalSed = Sed(wavelen=bulgeSed.wavelen, flambda=fl)

            sedList = [totalSed, bulgeSed, diskSed, agnSed]

            for i, spectrum in enumerate(sedList):
                if i==0:
                    msgroot = 'failed on total'
                elif i==1:
                    msgroot = 'failed on bulge'
                elif i==2:
                    msgroot = 'failed on disk'
                elif i==3:
                    msgroot = 'failed on agn'

                for j, b in enumerate(self.bandpasses):
                    controlSigma = calcMagError_sed(spectrum, self.totalBandpasses[j],
                                             self.skySeds[j],
                                             self.hardwareBandpasses[j],
                                             seeing=lsstDefaults.seeing(b),
                                             photParams=PhotometricParameters())

                    testSigma = line[26+(i*6)+j]
                    msg = '%e neq %e; ' % (testSigma, controlSigma) + msgroot
                    self.assertAlmostEqual(testSigma, controlSigma, 10, msg=msg)
                    ct += 1

        self.assertGreater(ct, 0)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(testPhotometricUncertaintyGetters)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
