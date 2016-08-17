from __future__ import with_statement
import numpy as np

import os
import unittest
import lsst.utils.tests
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.utils import myTestGals, myTestStars, \
                                     makeStarTestDB, makeGalTestDB, getOneChunk

from lsst.utils import getPackageDir
from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import PhotometricParameters, calcSNR_m5, LSSTdefaults
from lsst.sims.photUtils import calcSNR_sed, magErrorFromSNR
from lsst.sims.photUtils.utils import setM5
from lsst.sims.catUtils.utils import cartoonStars, cartoonGalaxies, testStars, testGalaxies, \
                                     cartoonStarsOnlyI, cartoonStarsIZ, \
                                     cartoonGalaxiesIG, galaxiesWithHoles
from lsst.sims.catUtils.mixins import PhotometryBase, PhotometryStars, PhotometryGalaxies


def setup_module(module):
    lsst.utils.tests.init()


class variabilityUnitTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create test databases
        if os.path.exists('PhotometryTestDatabase.db'):
            print "deleting database"
            os.unlink('PhotometryTestDatabase.db')

        makeStarTestDB(filename='PhotometryTestDatabase.db', size=100000, seedVal=1)
        makeGalTestDB(filename='PhotometryTestDatabase.db', size=100000, seedVal=1)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('PhotometryTestDatabase.db'):
            os.unlink('PhotometryTestDatabase.db')

    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.7,
                            boundType = 'circle',pointingRA=200.0,pointingDec=-30.0,
                            boundLength=1.0,
                            m5=[23.9, 25.0, 24.7, 24.0, 23.3, 22.1],
                            bandpassName=['u', 'g', 'r', 'i', 'z', 'y'])

        self.galaxy = myTestGals(database='PhotometryTestDatabase.db')
        self.star = myTestStars(database='PhotometryTestDatabase.db')

    def tearDown(self):
        del self.galaxy
        del self.star
        del self.obs_metadata

    def testGalaxyVariability(self):

        galcat = testGalaxies(self.galaxy, obs_metadata=self.obs_metadata)
        results = self.galaxy.query_columns(['varParamStr'], obs_metadata=self.obs_metadata,
                                             constraint='VarParamStr is not NULL')
        result = getOneChunk(results)
        ct = 0
        for row in result:
            mags=galcat.applyVariability(row['varParamStr'])
            ct += 1
        self.assertGreater(ct, 0) #to make sure that the test was actually performed

    def testStarVariability(self):
        starcat = testStars(self.star, obs_metadata=self.obs_metadata)
        results = self.star.query_columns(['varParamStr'], obs_metadata=self.obs_metadata,
                                         constraint='VarParamStr is not NULL')
        result = getOneChunk(results)
        ct = 0
        for row in result:
            ct += 1
            mags=starcat.applyVariability(row['varParamStr'])
        self.assertGreater(ct, 0) #to make sure that the test was actually performed

class photometryUnitTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create test databases
        if os.path.exists('PhotometryTestDatabase.db'):
            print "deleting database"
            os.unlink('PhotometryTestDatabase.db')

        makeStarTestDB(filename='PhotometryTestDatabase.db', size=100000, seedVal=1)
        makeGalTestDB(filename='PhotometryTestDatabase.db', size=100000, seedVal=1)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('PhotometryTestDatabase.db'):
            os.unlink('PhotometryTestDatabase.db')

    def setUp(self):
        defaults = LSSTdefaults()
        bandpassName=['u', 'g', 'r', 'i', 'z', 'y']
        self.obs_metadata = ObservationMetaData(mjd=52000.7,
                            bandpassName=bandpassName,
                            m5=[defaults.m5(mm) for mm in bandpassName],
                            boundType='circle',pointingRA=200.0,pointingDec=-30.0,
                            boundLength=1.0)

        self.galaxy = myTestGals(database='PhotometryTestDatabase.db')
        self.star = myTestStars(database='PhotometryTestDatabase.db')

    def tearDown(self):
        del self.galaxy
        del self.star
        del self.obs_metadata

    def testStarCatalog(self):
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests',
                               'scratchSpace', 'testPhotMixTestStarCat.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        test_cat=testStars(self.star, obs_metadata=self.obs_metadata)
        test_cat.write_catalog(catName)
        cat = open(catName)
        lines = cat.readlines()
        self.assertGreater(len(lines), 1) #to make sure we did not write an empty catalog
        cat.close()
        results = self.star.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)
        self.assertGreater(len(result), 0) #to make sure some results are returned
        if os.path.exists(catName):
            os.unlink(catName)


    def testGalaxyCatalog(self):
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests',
                               'scratchSpace', 'testPhotMixTestGalCat.txt')

        if os.path.exists(catName):
            os.unlink(catName)
        test_cat=testGalaxies(self.galaxy, obs_metadata=self.obs_metadata)
        test_cat.write_catalog(catName)
        cat = open(catName)
        lines = cat.readlines()
        self.assertGreater(len(lines), 1) #to make sure we did not write an empty catalog
        cat.close()
        results = self.galaxy.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)
        self.assertGreater(len(result), 0) #to make sure some results are returned

        if os.path.exists(catName):
            os.unlink(catName)


    def test_m5_exceptions(self):
        """
        Test that the correct exception is raised when you ask for a photometric
        uncertainty but do not define the required m5 value
        """

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'm5_exception_cat.txt')

        obs = ObservationMetaData(pointingRA=25.0, pointingDec=-14.0,
                                  boundType='circle', boundLength=0.1,
                                  bandpassName=['u', 'g', 'r', 'z', 'y'],
                                  m5 = [24.0] * 5,
                                  mjd=57388.0)

        with self.assertRaises(KeyError) as context:
            cat = testStars(self.star, obs_metadata=obs)
            cat.write_catalog(catName)

        self.assertIn('Is it possible your ObservationMetaData does not have the proper\nm5 values defined?',
                      context.exception.args[0])


        with self.assertRaises(KeyError) as context:
            cat = testGalaxies(self.galaxy, obs_metadata=obs)
            cat.write_catalog(catName)

        self.assertIn('Is it possible your ObservationMetaData does not have the proper\nm5 values defined?',
                      context.exception.args[0])

        if os.path.exists(catName):
            os.unlink(catName)


    def testSumMagnitudes(self):
        """
        Test that the method sum_magnitudes in PhotometryGalaxies handles
        NaNs correctly.  Test it both in the vectorized and non-vectorized form.
        """
        mm_0 = 22.0

        bulge = 15.0*np.ones(8)

        disk = 15.2*np.ones(8)

        agn = 15.4*np.ones(8)

        bulge[0] = np.NaN
        disk[1] = np.NaN
        agn[2] = np.NaN

        bulge[3] = np.NaN
        disk[3] = np.NaN

        bulge[4] = np.NaN
        agn[4] = np.NaN

        disk[5] = np.NaN
        agn[5] = np.NaN

        bulge[7] = np.NaN
        disk[7] = np.NaN
        agn[7] = np.NaN

        bulge_flux = np.power(10.0, -0.4*(bulge-mm_0))
        disk_flux = np.power(10.0, -0.4*(disk-mm_0))
        agn_flux = np.power(10.0, -0.4*(agn-mm_0))

        answer = np.zeros(8)
        answer[0] = -2.5*np.log10(disk_flux[0]+agn_flux[0]) + mm_0
        answer[1] = -2.5*np.log10(bulge_flux[1]+agn_flux[1]) + mm_0
        answer[2] = -2.5*np.log10(bulge_flux[2]+disk_flux[2]) + mm_0
        answer[3] = -2.5*np.log10(agn_flux[3]) + mm_0
        answer[4] = -2.5*np.log10(disk_flux[4]) + mm_0
        answer[5] = -2.5*np.log10(bulge_flux[5]) + mm_0
        answer[6] = -2.5*np.log10(bulge_flux[6]+disk_flux[6]+agn_flux[6]) + mm_0
        answer[7] = np.NaN

        phot = PhotometryGalaxies()
        test = phot.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

        np.testing.assert_array_almost_equal(test, answer, decimal=10)

        for ix, (bb, dd, aa, truth) in enumerate(zip(bulge, disk, agn, answer)):
            test = phot.sum_magnitudes(bulge=bb, disk=dd, agn=aa)
            if ix<7:
                self.assertAlmostEqual(test, truth, 10)
                self.assertFalse(np.isnan(test))
            else:
                self.assertTrue(np.isnan(test))
                self.assertTrue(np.isnan(truth))

    def testSumMagnitudesCatalog(self):
        """
        test that sum_magnitudes handles NaNs correctly in the context
        of a catalog by outputting a catalog of galaxies with NaNs in
        different component magnitudes, reading that catalog back in,
        and then calculating the summed magnitude by hand and comparing
        """

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'testPhotMixGalaxiesWithHoles.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        obs_metadata=ObservationMetaData(mjd=50000.0,
                               boundType='circle',pointingRA=0.0,pointingDec=0.0,
                               boundLength=10.0)
        test_cat=galaxiesWithHoles(self.galaxy,obs_metadata=obs_metadata)
        test_cat.write_catalog(catName)

        dtype = np.dtype([
                            ('raJ2000', np.float),
                            ('decJ2000', np.float),
                            ('u', np.float), ('g', np.float), ('r', np.float),
                            ('i', np.float), ('z', np.float), ('y', np.float),
                            ('ub', np.float), ('gb', np.float), ('rb', np.float),
                            ('ib', np.float), ('zb', np.float), ('yb', np.float),
                            ('ud', np.float), ('gd', np.float), ('rd', np.float),
                            ('id', np.float), ('zd', np.float), ('yd', np.float),
                            ('ua', np.float), ('ga', np.float), ('ra', np.float),
                            ('ia', np.float), ('za', np.float), ('ya', np.float)
                            ])



        data = np.genfromtxt(catName, dtype=dtype, delimiter=', ')
        self.assertGreater(len(data), 16)
        phot = PhotometryGalaxies()

        test = phot.sum_magnitudes(bulge=data['ub'], disk=data['ud'], agn=data['ua'])
        np.testing.assert_array_almost_equal(test, data['u'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['gb'], disk=data['gd'], agn=data['ga'])
        np.testing.assert_array_almost_equal(test, data['g'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['rb'], disk=data['rd'], agn=data['ra'])
        np.testing.assert_array_almost_equal(test, data['r'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['ib'], disk=data['id'], agn=data['ia'])
        np.testing.assert_array_almost_equal(test, data['i'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['zb'], disk=data['zd'], agn=data['za'])
        np.testing.assert_array_almost_equal(test, data['z'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['yb'], disk=data['yd'], agn=data['ya'])
        np.testing.assert_array_almost_equal(test, data['y'], decimal=10)

        # make sure that there were some NaNs for our catalog to deal with (but that they were not
        # all NaNs
        for line in [data['u'], data['g'], data['r'], data['i'], data['z'], data['y'],
                     data['ub'], data['gb'], data['rb'], data['ib'], data['zb'], data['yb'],
                     data['ud'], data['gd'], data['rd'], data['id'], data['zd'], data['yd'],
                     data['ua'], data['ga'], data['ra'], data['ia'], data['za'], data['ya']]:

            ctNans = len(np.where(np.isnan(line))[0])
            self.assertGreater(ctNans, 0)
            self.assertLess(ctNans, len(line))

        if os.path.exists(catName):
            os.unlink(catName)


    def testAlternateBandpassesStars(self):
        """
        This will test our ability to do photometry using non-LSST bandpasses.

        It will first calculate the magnitudes using the getters in cartoonPhotometryStars.

        It will then load the alternate bandpass files 'by hand' and re-calculate the magnitudes
        and make sure that the magnitude values agree.  This is guarding against the possibility
        that some default value did not change and the code actually ended up loading the
        LSST bandpasses.
        """

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'testPhotMixTestStarsCartoon.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        obs_metadata_pointed=ObservationMetaData(mjd=2013.23,
                                                 boundType='circle',pointingRA=200.0,pointingDec=-30.0,
                                                 boundLength=1.0)

        test_cat=cartoonStars(self.star,obs_metadata=obs_metadata_pointed)

        test_cat.write_catalog(catName)

        cartoonDir = os.path.join(getPackageDir('sims_photUtils'), 'tests', 'cartoonSedTestData')
        testBandPasses = {}
        keys = ['u','g','r','i','z']

        bplist = []

        for kk in keys:
            testBandPasses[kk] = Bandpass()
            testBandPasses[kk].readThroughput(os.path.join(cartoonDir,"test_bandpass_%s.dat" % kk))
            bplist.append(testBandPasses[kk])

        sedObj = Sed()
        phiArray, waveLenStep = sedObj.setupPhiArray(bplist)

        i = 0

        #since all of the SEDs in the cartoon database are the same, just test on the first
        #if we ever include more SEDs, this can be something like
        #for ss in test_cata.sedMasterList:
        #
        ss=test_cat.sedMasterList[0]
        ss.resampleSED(wavelen_match = bplist[0].wavelen)
        ss.flambdaTofnu()
        mags = -2.5*np.log10(np.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
        self.assertEqual(len(mags), len(test_cat.cartoonBandpassDict))
        self.assertGreater(len(mags), 0)
        for j in range(len(mags)):
            self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterList[i][j],4)

        with open(catName, 'r') as input_file:
            lines = input_file.readlines()
            self.assertGreater(len(lines), 1)

        if os.path.exists(catName):
            os.unlink(catName)

    def testAlternateBandpassesGalaxies(self):
        """
        the same as testAlternateBandpassesStars, but for galaxies
        """

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'testPhotMix_testAlternateBandpassesGalaxies.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        obs_metadata_pointed=ObservationMetaData(mjd=50000.0,
                               boundType='circle',pointingRA=0.0,pointingDec=0.0,
                               boundLength=10.0)

        test_cat=cartoonGalaxies(self.galaxy,obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(catName)

        dtype = np.dtype([
                             ('galid', np.int),
                             ('ra', np.float),
                             ('dec', np.float),
                             ('uTotal', np.float),
                             ('gTotal', np.float),
                             ('rTotal', np.float),
                             ('iTotal', np.float),
                             ('zTotal', np.float),
                             ('uBulge', np.float),
                             ('gBulge', np.float),
                             ('rBulge', np.float),
                             ('iBulge', np.float),
                             ('zBulge', np.float),
                             ('uDisk', np.float),
                             ('gDisk', np.float),
                             ('rDisk', np.float),
                             ('iDisk', np.float),
                             ('zDisk', np.float),
                             ('uAgn', np.float),
                             ('gAgn', np.float),
                             ('rAgn', np.float),
                             ('iAgn', np.float),
                             ('zAgn', np.float),
                             ('bulgeName', str, 200),
                             ('bulgeNorm', np.float),
                             ('bulgeAv', np.float),
                             ('diskName', str, 200),
                             ('diskNorm', np.float),
                             ('diskAv', np.float),
                             ('agnName', str, 200),
                             ('agnNorm', np.float),
                             ('redshift', np.float)
                            ])

        catData = np.genfromtxt(catName, dtype=dtype, delimiter=', ')

        self.assertGreater(len(catData), 0)

        cartoonDir = getPackageDir('sims_photUtils')
        cartoonDir = os.path.join(cartoonDir, 'tests', 'cartoonSedTestData')
        sedDir = getPackageDir('sims_sed_library')

        testBandpasses = {}
        keys = ['u','g','r','i','z']

        for kk in keys:
            testBandpasses[kk] = Bandpass()
            testBandpasses[kk].readThroughput(os.path.join(cartoonDir,"test_bandpass_%s.dat" % kk))

        imsimBand = Bandpass()
        imsimBand.imsimBandpass()

        specMap = defaultSpecMap

        ct = 0
        for line in catData:
            bulgeMagList = []
            diskMagList = []
            agnMagList = []
            if line['bulgeName'] == 'None':
                for bp in keys:
                    self.assertTrue(np.isnan(line['%sBulge' % bp]))
                    bulgeMagList.append(np.NaN)
            else:
                ct += 1
                dummySed = Sed()
                dummySed.readSED_flambda(os.path.join(sedDir, specMap[line['bulgeName']]))
                fnorm = dummySed.calcFluxNorm(line['bulgeNorm'], imsimBand)
                dummySed.multiplyFluxNorm(fnorm)
                a_int, b_int = dummySed.setupCCMab()
                dummySed.addCCMDust(a_int, b_int, A_v=line['bulgeAv'])
                dummySed.redshiftSED(line['redshift'], dimming=True)
                dummySed.resampleSED(wavelen_match=testBandpasses['u'].wavelen)
                for bpName in keys:
                    mag = dummySed.calcMag(testBandpasses[bpName])
                    self.assertAlmostEqual(mag, line['%sBulge' % bpName], 10)
                    bulgeMagList.append(mag)

            if line['diskName'] == 'None':
                for bp in keys:
                    self.assertTrue(np.isnan(line['%sDisk' % bp]))
                    diskMagList.append(np.NaN)
            else:
                ct += 1
                dummySed = Sed()
                dummySed.readSED_flambda(os.path.join(sedDir, specMap[line['diskName']]))
                fnorm = dummySed.calcFluxNorm(line['diskNorm'], imsimBand)
                dummySed.multiplyFluxNorm(fnorm)
                a_int, b_int = dummySed.setupCCMab()
                dummySed.addCCMDust(a_int, b_int, A_v=line['diskAv'])
                dummySed.redshiftSED(line['redshift'], dimming=True)
                dummySed.resampleSED(wavelen_match=testBandpasses['u'].wavelen)
                for bpName in keys:
                    mag = dummySed.calcMag(testBandpasses[bpName])
                    self.assertAlmostEqual(mag, line['%sDisk' % bpName], 10)
                    diskMagList.append(mag)

            if line['agnName'] == 'None':
                for bp in keys:
                    self.assertTrue(np.isnan(line['%sAgn' % bp]))
                    agnMagList.append(np.NaN)
            else:
                ct += 1
                dummySed = Sed()
                dummySed.readSED_flambda(os.path.join(sedDir, specMap[line['agnName']]))
                fnorm = dummySed.calcFluxNorm(line['agnNorm'], imsimBand)
                dummySed.multiplyFluxNorm(fnorm)
                dummySed.redshiftSED(line['redshift'], dimming=True)
                dummySed.resampleSED(wavelen_match=testBandpasses['u'].wavelen)
                for bpName in keys:
                    mag = dummySed.calcMag(testBandpasses[bpName])
                    self.assertAlmostEqual(mag, line['%sAgn' % bpName], 10)
                    agnMagList.append(mag)

            totalMags = PhotometryGalaxies().sum_magnitudes(bulge=np.array(bulgeMagList),
                                                            disk=np.array(diskMagList),
                                                            agn=np.array(agnMagList))

            for testMag, bpName in zip(totalMags, keys):
                if np.isnan(line['%sTotal' % bpName]):
                    self.assertTrue(np.isnan(testMag))
                else:
                    self.assertAlmostEqual(testMag, line['%sTotal' % bpName],10)


        self.assertGreater(ct, 0)
        if os.path.exists(catName):
            os.unlink(catName)

    def testStellarPhotometryIndices(self):
        """
        A test to make sure that stellar photometry still calculates the right values
        even when it is not calculating all of the magnitudes in the getter
        """

        baselineDtype = np.dtype([('id',int),
                                     ('raObserved', float), ('decObserved', float),
                                     ('magNorm', float),
                                     ('cartoon_u', float), ('cartoon_g',float),
                                     ('cartoon_r', float), ('cartoon_i', float),
                                     ('cartoon_z', float)])

        baselineCatName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                      'testPhotMix_testStellarIndices_stellarBaselineCatalog.txt')

        if os.path.exists(baselineCatName):
            os.unlink(baselineCatName)

        testDtype = np.dtype([('id',int),
                                 ('raObserved',float), ('decObserved',float),
                                 ('cartoon_i',float)])

        testCatName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                  'testPhotMix_testStellarIndices_stellarTestCatalog.txt')

        if os.path.exists(testCatName):
            os.unlink(testCatName)

        obs_metadata_pointed=ObservationMetaData(mjd=2013.23,
                                                 boundType='circle',pointingRA=200.0,pointingDec=-30.0,
                                                 boundLength=1.0)

        baseline_cat=cartoonStars(self.star,obs_metadata=obs_metadata_pointed)
        baseline_cat.write_catalog(baselineCatName)
        baselineData = np.genfromtxt(baselineCatName, dtype=baselineDtype, delimiter=',')
        self.assertGreater(len(baselineData), 0)

        test_cat=cartoonStarsOnlyI(self.star, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = np.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        self.assertGreater(len(testData), 0)

        for b, t in zip(baselineData, testData):
            self.assertAlmostEqual(b['cartoon_i'], t['cartoon_i'], 10)

        testDtype = np.dtype([('id',int),
                                 ('raObserved',float), ('decObserved',float),
                                 ('cartoon_i',float), ('cartoon_z',float)])


        test_cat=cartoonStarsIZ(self.star, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = np.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        self.assertGreater(len(testData), 0)

        for b, t in zip(baselineData, testData):
            self.assertAlmostEqual(b['cartoon_i'], t['cartoon_i'], 10)
            self.assertAlmostEqual(b['cartoon_z'], t['cartoon_z'], 10)

        if os.path.exists(testCatName):
            os.unlink(testCatName)
        if os.path.exists(baselineCatName):
            os.unlink(baselineCatName)

    def testGalaxyPhotometricIndices(self):
        baselineCatName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                       'testPhotMix_testGalIndices_galaxyBaselineCatalog.txt')

        if os.path.exists(baselineCatName):
            os.unlink(baselineCatName)

        baselineDtype = np.dtype([('galid', int),
                                     ('raObserved', float),
                                     ('decObserved', float),
                                     ('ctotal_u', float),
                                     ('ctotal_g', float),
                                     ('ctotal_r', float),
                                     ('ctotal_i', float),
                                     ('ctotal_z', float)])

        obs_metadata_pointed=ObservationMetaData(mjd=50000.0,
                               boundType='circle',pointingRA=0.0,pointingDec=0.0,
                               boundLength=10.0)

        baseline_cat=cartoonGalaxies(self.galaxy,obs_metadata=obs_metadata_pointed)
        baseline_cat.write_catalog(baselineCatName)
        baselineData = np.genfromtxt(baselineCatName, dtype=baselineDtype, delimiter=',')
        self.assertGreater(len(baselineData), 0)

        testCatName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                   'testPhotMix_testGalIndices_galaxyTestCatalog.txt')

        if os.path.exists(testCatName):
            os.unlink(testCatName)

        testDtype = np.dtype([('galid', int),
                                 ('raObserved', float),
                                 ('decObserved', float),
                                 ('ctotal_i', float),
                                 ('ctotal_g', float)])

        test_cat = cartoonGalaxiesIG(self.galaxy, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = np.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        self.assertGreater(len(testData), 0)

        for b,t in zip(baselineData, testData):
            self.assertAlmostEqual(b['ctotal_i'], t['ctotal_i'], 10)
            self.assertAlmostEqual(b['ctotal_g'], t['ctotal_g'], 10)

        if os.path.exists(baselineCatName):
            os.unlink(baselineCatName)

        if os.path.exists(testCatName):
            os.unlink(testCatName)

    def testPhotometricIndicesRaw(self):
        """
        Use manMagCalc_list with specified indices on an Sed.  Make sure
        that the appropriate magnitudes are or are not Nan
        """
        starName = os.path.join(getPackageDir('sims_sed_library'),defaultSpecMap['km20_5750.fits_g40_5790'])
        starPhot = BandpassDict.loadTotalBandpassesFromFiles()
        testSed = Sed()
        testSed.readSED_flambda(starName)
        indices = [1,3]
        mags = starPhot.magListForSed(testSed, indices=indices)
        self.assertTrue(np.isnan(mags[0]))
        self.assertFalse(np.isnan(mags[1]))
        self.assertTrue(np.isnan(mags[2]))
        self.assertFalse(np.isnan(mags[3]))
        self.assertTrue(np.isnan(mags[4]))
        self.assertTrue(np.isnan(mags[5]))
        self.assertEqual(len(mags), 6)


class UncertaintyMixinTest(unittest.TestCase):
    def setUp(self):
        starName = os.path.join(getPackageDir('sims_sed_library'),defaultSpecMap['km20_5750.fits_g40_5790'])
        self.starSED = Sed()
        self.starSED.readSED_flambda(starName)
        imsimband = Bandpass()
        imsimband.imsimBandpass()
        fNorm = self.starSED.calcFluxNorm(22.0, imsimband)
        self.starSED.multiplyFluxNorm(fNorm)

        self.totalBandpasses = []
        self.hardwareBandpasses = []

        componentList = ['detector.dat', 'm1.dat', 'm2.dat', 'm3.dat',
                         'lens1.dat', 'lens2.dat', 'lens3.dat']
        hardwareComponents = []
        for c in componentList:
            hardwareComponents.append(os.path.join(getPackageDir('throughputs'),'baseline',c))

        self.bandpasses = ['u', 'g', 'r', 'i', 'z', 'y']
        for b in self.bandpasses:
            filterName = os.path.join(getPackageDir('throughputs'),'baseline','filter_%s.dat' % b)
            components = hardwareComponents + [filterName]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.hardwareBandpasses.append(bandpassDummy)
            components = components + [os.path.join(getPackageDir('throughputs'),'baseline','atmos.dat')]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.totalBandpasses.append(bandpassDummy)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
