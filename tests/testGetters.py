import os
import numpy
import unittest
import eups

import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB
from lsst.sims.catalogs.measures.instance import InstanceCatalog, defaultSpecMap
from lsst.sims.catUtils.utils import testStarsDBObj, testGalaxyTileDBObj
from lsst.sims.photUtils import Sed, Bandpass, PhotometricDefaults, setM5
from lsst.sims.photUtils import PhotometryStars, PhotometryGalaxies

class testStarCatalog(InstanceCatalog, PhotometryStars):
    sig2sys = 0.0003

    column_outputs = ['raJ2000', 'decJ2000',
                      'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r', 'sigma_lsst_i',
                      'sigma_lsst_z', 'sigma_lsst_y', 'sedFilename', 'magNorm']

class testGalaxyCatalog(InstanceCatalog, PhotometryGalaxies):
    sig2sys = 0.0003

    column_outputs = ['raJ2000', 'decJ2000',
                      'uRecalc', 'gRecalc', 'rRecalc', 'iRecalc', 'zRecalc', 'yRecalc',
                      'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
                      'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
                      'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn',
                      'sigma_uRecalc', 'sigma_gRecalc', 'sigma_rRecalc', 'sigma_iRecalc',
                      'sigma_zRecalc', 'sigma_yRecalc',
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
        cls.dbName = 'uncertaintyTestDB.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        cls.obs_metadata = makePhoSimTestDB(filename=cls.dbName, size=10, radius = 5.0)
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
            normalizedSedDummy = setM5(cls.obs_metadata.m5(cls.bandpasses[i]), sedDummy,
                                       cls.totalBandpasses[i], cls.hardwareBandpasses[i],
                                       seeing=PhotometricDefaults.seeing[cls.bandpasses[i]])
            cls.skySeds.append(normalizedSedDummy)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)
        del cls.dbName
        del cls.connectionString
        del cls.obs_metadata
        del cls.totalBandpasses
        del cls.hardwareBandpasses
        del cls.skySeds

    def testStellarPhotometricUncertainties(self):
        """
        Test in the case of a catalog of stars
        """
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

    def testGalaxyPhotometricUncertainties(self):
        """
        Test in the case of a catalog of galaxies
        """
        phot = PhotometryGalaxies()
        phot.loadTotalBandpassesFromFiles()
        galDB = testGalaxyTileDBObj(address=self.connectionString)
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
            bulgeSed.readSED_flambda(os.path.join(eups.productDir('sims_sed_library'),
                                     defaultSpecMap[bulgeSedName]))
            fNorm=bulgeSed.calcFluxNorm(magNormBulge, imsimband)
            bulgeSed.multiplyFluxNorm(fNorm)

            diskSed = Sed()
            diskSed.readSED_flambda(os.path.join(eups.productDir('sims_sed_library'),
                                    defaultSpecMap[diskSedName]))
            fNorm = diskSed.calcFluxNorm(magNormDisk, imsimband)
            diskSed.multiplyFluxNorm(fNorm)

            agnSed = Sed()
            agnSed.readSED_flambda(os.path.join(eups.productDir('sims_sed_library'),
                                   defaultSpecMap[agnSedName]))
            fNorm = agnSed.calcFluxNorm(magNormAgn, imsimband)
            agnSed.multiplyFluxNorm(fNorm)


            phot.applyAvAndRedshift([bulgeSed, diskSed], internalAv = [avBulge, avDisk], redshift=[redshift, redshift])
            phot.applyAvAndRedshift([agnSed], redshift=[redshift])

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
                    controlSNR = spectrum.calcSNR_psf(self.totalBandpasses[j],
                                                      self.skySeds[j],
                                                      self.hardwareBandpasses[j],
                                                      seeing=PhotometricDefaults.seeing[b])

                    controlNoverS = 1./(controlSNR*controlSNR) + galCat.sig2sys
                    controlSigma = 2.5*numpy.log10(1.0+numpy.sqrt(controlNoverS))
                    testSigma = line[26+(i*6)+j]
                    msg = '%e neq %e; ' % (testSigma, controlSigma) + msgroot
                    self.assertAlmostEqual(testSigma, controlSigma, 10, msg=msg)
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
