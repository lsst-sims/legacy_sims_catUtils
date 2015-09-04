import numpy

import os
import unittest
import lsst.utils
import lsst.utils.tests as utilsTests
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.generation.utils import myTestGals, myTestStars, \
                                                makeStarTestDB, makeGalTestDB, getOneChunk

from lsst.sims.utils import defaultSpecMap
from lsst.sims.photUtils.Bandpass import Bandpass
from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils import loadTotalBandpassesFromFiles
from lsst.sims.catUtils.utils import cartoonStars, cartoonGalaxies, testStars, testGalaxies, \
                                     cartoonStarsOnlyI, cartoonStarsIZ, \
                                     cartoonGalaxiesIG, galaxiesWithHoles
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies


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
                            boundType = 'circle',unrefractedRA=200.0,unrefractedDec=-30.0,
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
        self.assertTrue(ct>0) #to make sure that the test was actually performed

    def testStarVariability(self):
        starcat = testStars(self.star, obs_metadata=self.obs_metadata)
        results = self.star.query_columns(['varParamStr'], obs_metadata=self.obs_metadata,
                                         constraint='VarParamStr is not NULL')
        result = getOneChunk(results)
        ct = 0
        for row in result:
            ct += 1
            mags=starcat.applyVariability(row['varParamStr'])
        self.assertTrue(ct>0) #to make sure that the test was actually performed

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
        self.obs_metadata = ObservationMetaData(mjd=52000.7, bandpassName='i',
                            boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                            boundLength=1.0, m5 = 25.0)

        self.galaxy = myTestGals(database='PhotometryTestDatabase.db')
        self.star = myTestStars(database='PhotometryTestDatabase.db')

    def tearDown(self):
        del self.galaxy
        del self.star
        del self.obs_metadata

    def testStarCatalog(self):
        test_cat=testStars(self.star, obs_metadata=self.obs_metadata)
        test_cat.write_catalog("testStarsOutput.txt")
        cat = open("testStarsOutput.txt")
        lines = cat.readlines()
        self.assertTrue(len(lines)>1) #to make sure we did not write an empty catalog
        cat.close()
        results = self.star.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)
        self.assertTrue(len(result)>0) #to make sure some results are returned
        os.unlink("testStarsOutput.txt")


    def testGalaxyCatalog(self):
        test_cat=testGalaxies(self.galaxy, obs_metadata=self.obs_metadata)
        test_cat.write_catalog("testGalaxiesOutput.txt")
        cat = open("testGalaxiesOutput.txt")
        lines = cat.readlines()
        self.assertTrue(len(lines)>1) #to make sure we did not write an empty catalog
        cat.close()
        results = self.galaxy.query_columns(obs_metadata=self.obs_metadata)
        result = getOneChunk(results)
        self.assertTrue(len(result)>0) #to make sure some results are returned
        os.unlink("testGalaxiesOutput.txt")

    def testSumMagnitudes(self):
        """
        Test that the method sum_magnitudes in PhotometryGalaxies handles
        NaNs correctly.  Test it both in the vectorized and non-vectorized form.
        """
        mm_0 = 22.0

        bulge = 15.0*numpy.ones(8)

        disk = 15.2*numpy.ones(8)

        agn = 15.4*numpy.ones(8)

        bulge[0] = numpy.NaN
        disk[1] = numpy.NaN
        agn[2] = numpy.NaN

        bulge[3] = numpy.NaN
        disk[3] = numpy.NaN

        bulge[4] = numpy.NaN
        agn[4] = numpy.NaN

        disk[5] = numpy.NaN
        agn[5] = numpy.NaN

        bulge[7] = numpy.NaN
        disk[7] = numpy.NaN
        agn[7] = numpy.NaN

        bulge_flux = numpy.power(10.0, -0.4*(bulge-mm_0))
        disk_flux = numpy.power(10.0, -0.4*(disk-mm_0))
        agn_flux = numpy.power(10.0, -0.4*(agn-mm_0))

        answer = numpy.zeros(8)
        answer[0] = -2.5*numpy.log10(disk_flux[0]+agn_flux[0]) + mm_0
        answer[1] = -2.5*numpy.log10(bulge_flux[1]+agn_flux[1]) + mm_0
        answer[2] = -2.5*numpy.log10(bulge_flux[2]+disk_flux[2]) + mm_0
        answer[3] = -2.5*numpy.log10(agn_flux[3]) + mm_0
        answer[4] = -2.5*numpy.log10(disk_flux[4]) + mm_0
        answer[5] = -2.5*numpy.log10(bulge_flux[5]) + mm_0
        answer[6] = -2.5*numpy.log10(bulge_flux[6]+disk_flux[6]+agn_flux[6]) + mm_0
        answer[7] = numpy.NaN

        phot = PhotometryGalaxies()
        test = phot.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

        numpy.testing.assert_array_almost_equal(test, answer, decimal=10)

        for ix, (bb, dd, aa, truth) in enumerate(zip(bulge, disk, agn, answer)):
            test = phot.sum_magnitudes(bulge=bb, disk=dd, agn=aa)
            if ix<7:
                self.assertAlmostEqual(test, truth, 10)
                self.assertTrue(not numpy.isnan(test))
            else:
                self.assertTrue(numpy.isnan(test))
                self.assertTrue(numpy.isnan(truth))

    def testSumMagnitudesCatalog(self):
        """
        test that sum_magnitudes handles NaNs correctly in the context
        of a catalog by outputting a catalog of galaxies with NaNs in
        different component magnitudes, reading that catalog back in,
        and then calculating the summed magnitude by hand and comparing
        """

        catName = 'galaxiesWithHoles.txt'
        obs_metadata=ObservationMetaData(mjd=50000.0,
                               boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                               boundLength=10.0)
        test_cat=galaxiesWithHoles(self.galaxy,obs_metadata=obs_metadata)
        test_cat.write_catalog(catName)

        dtype = numpy.dtype([
                            ('raJ2000', numpy.float),
                            ('decJ2000', numpy.float),
                            ('u', numpy.float), ('g', numpy.float), ('r', numpy.float),
                            ('i', numpy.float), ('z', numpy.float), ('y', numpy.float),
                            ('ub', numpy.float), ('gb', numpy.float), ('rb', numpy.float),
                            ('ib', numpy.float), ('zb', numpy.float), ('yb', numpy.float),
                            ('ud', numpy.float), ('gd', numpy.float), ('rd', numpy.float),
                            ('id', numpy.float), ('zd', numpy.float), ('yd', numpy.float),
                            ('ua', numpy.float), ('ga', numpy.float), ('ra', numpy.float),
                            ('ia', numpy.float), ('za', numpy.float), ('ya', numpy.float)
                            ])



        data = numpy.genfromtxt(catName, dtype=dtype, delimiter=', ')
        self.assertTrue(len(data)>16)
        phot = PhotometryGalaxies()

        test = phot.sum_magnitudes(bulge=data['ub'], disk=data['ud'], agn=data['ua'])
        numpy.testing.assert_array_almost_equal(test, data['u'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['gb'], disk=data['gd'], agn=data['ga'])
        numpy.testing.assert_array_almost_equal(test, data['g'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['rb'], disk=data['rd'], agn=data['ra'])
        numpy.testing.assert_array_almost_equal(test, data['r'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['ib'], disk=data['id'], agn=data['ia'])
        numpy.testing.assert_array_almost_equal(test, data['i'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['zb'], disk=data['zd'], agn=data['za'])
        numpy.testing.assert_array_almost_equal(test, data['z'], decimal=10)

        test = phot.sum_magnitudes(bulge=data['yb'], disk=data['yd'], agn=data['ya'])
        numpy.testing.assert_array_almost_equal(test, data['y'], decimal=10)

        # make sure that there were some NaNs for our catalog to deal with (but that they were not
        # all NaNs
        for line in [data['u'], data['g'], data['r'], data['i'], data['z'], data['y'],
                     data['ub'], data['gb'], data['rb'], data['ib'], data['zb'], data['yb'],
                     data['ud'], data['gd'], data['rd'], data['id'], data['zd'], data['yd'],
                     data['ua'], data['ga'], data['ra'], data['ia'], data['za'], data['ya']]:

            ctNans = len(numpy.where(numpy.isnan(line))[0])
            self.assertTrue(ctNans>0)
            self.assertTrue(ctNans<len(line))

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

        obs_metadata_pointed=ObservationMetaData(mjd=2013.23,
                                                 boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                                                 boundLength=1.0)

        test_cat=cartoonStars(self.star,obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog("testStarsCartoon.txt")

        cartoonDir = os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'
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
        mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
        self.assertTrue(len(mags)==len(test_cat.bandpassDict))
        self.assertTrue(len(mags)>0)
        for j in range(len(mags)):
            self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterList[i][j],10)
        i += 1

        os.unlink("testStarsCartoon.txt")

    def testAlternateBandpassesGalaxies(self):
        """
        the same as testAlternateBandpassesStars, but for galaxies
        """

        obs_metadata_pointed=ObservationMetaData(mjd=50000.0,
                               boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                               boundLength=10.0)

        test_cat=cartoonGalaxies(self.galaxy,obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog("testGalaxiesCartoon.txt")

        cartoonDir = os.getenv('SIMS_PHOTUTILS_DIR')+'/tests/cartoonSedTestData/'
        testBandPasses = {}
        keys = ['u','g','r','i','z']

        bplist = []

        for kk in keys:
            testBandPasses[kk] = Bandpass()
            testBandPasses[kk].readThroughput(os.path.join(cartoonDir,"test_bandpass_%s.dat" % kk))
            bplist.append(testBandPasses[kk])

        sedObj = Sed()
        phiArray, waveLenStep = sedObj.setupPhiArray(bplist)

        components = ['Bulge', 'Disk', 'Agn']

        ct = 0
        for cc in components:
            i = 0

            for ss in test_cat.sedMasterDict[cc]:
                if ss.wavelen != None:
                    ss.resampleSED(wavelen_match = bplist[0].wavelen)
                    ss.flambdaTofnu()
                    mags = -2.5*numpy.log10(numpy.sum(phiArray*ss.fnu, axis=1)*waveLenStep) - ss.zp
                    for j in range(len(mags)):
                        ct += 1
                        self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterDict[cc][i][j],10)
                i += 1

        self.assertTrue(ct>0)
        os.unlink("testGalaxiesCartoon.txt")

    def testStellarPhotometryIndices(self):
        """
        A test to make sure that stellar photometry still calculates the right values
        even when it is not calculating all of the magnitudes in the getter
        """

        baselineDtype = numpy.dtype([('id',int),
                                     ('raObserved', float), ('decObserved', float),
                                     ('magNorm', float),
                                     ('cartoon_u', float), ('cartoon_g',float),
                                     ('cartoon_r', float), ('cartoon_i', float),
                                     ('cartoon_z', float)])

        baselineCatName = 'stellarBaselineCatalog.txt'

        testDtype = numpy.dtype([('id',int),
                                 ('raObserved',float), ('decObserved',float),
                                 ('cartoon_i',float)])

        testCatName = 'stellarTestCatalog.txt'


        obs_metadata_pointed=ObservationMetaData(mjd=2013.23,
                                                 boundType='circle',unrefractedRA=200.0,unrefractedDec=-30.0,
                                                 boundLength=1.0)

        baseline_cat=cartoonStars(self.star,obs_metadata=obs_metadata_pointed)
        baseline_cat.write_catalog(baselineCatName)
        baselineData = numpy.genfromtxt(baselineCatName, dtype=baselineDtype, delimiter=',')

        test_cat=cartoonStarsOnlyI(self.star, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = numpy.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        ct = 0
        for b, t in zip(baselineData, testData):
            self.assertAlmostEqual(b['cartoon_i'], t['cartoon_i'], 10)
            ct+=1
        self.assertTrue(ct>0)

        testDtype = numpy.dtype([('id',int),
                                 ('raObserved',float), ('decObserved',float),
                                 ('cartoon_i',float), ('cartoon_z',float)])


        test_cat=cartoonStarsIZ(self.star, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = numpy.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        ct = 0
        for b, t in zip(baselineData, testData):
            self.assertAlmostEqual(b['cartoon_i'], t['cartoon_i'], 10)
            self.assertAlmostEqual(b['cartoon_z'], t['cartoon_z'], 10)
            ct+=1
        self.assertTrue(ct>0)

        if os.path.exists(testCatName):
            os.unlink(testCatName)
        if os.path.exists(baselineCatName):
            os.unlink(baselineCatName)

    def testGalaxyPhotometricIndices(self):
        baselineCatName = 'galaxyBaselineCatalog.txt'
        baselineDtype = numpy.dtype([('galid', int),
                                     ('raObserved', float),
                                     ('decObserved', float),
                                     ('ctotal_u', float),
                                     ('ctotal_g', float),
                                     ('ctotal_r', float),
                                     ('ctotal_i', float),
                                     ('ctotal_z', float)])

        obs_metadata_pointed=ObservationMetaData(mjd=50000.0,
                               boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                               boundLength=10.0)

        baseline_cat=cartoonGalaxies(self.galaxy,obs_metadata=obs_metadata_pointed)
        baseline_cat.write_catalog(baselineCatName)
        baselineData = numpy.genfromtxt(baselineCatName, dtype=baselineDtype, delimiter=',')

        testCatName = 'galaxyTestCatalog.txt'
        testDtype = numpy.dtype([('galid', int),
                                 ('raObserved', float),
                                 ('decObserved', float),
                                 ('ctotal_i', float),
                                 ('ctotal_g', float)])
        test_cat = cartoonGalaxiesIG(self.galaxy, obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(testCatName)
        testData = numpy.genfromtxt(testCatName, dtype=testDtype, delimiter=',')
        ct = 0
        for b,t in zip(baselineData, testData):
            self.assertAlmostEqual(b['ctotal_i'], t['ctotal_i'], 10)
            self.assertAlmostEqual(b['ctotal_g'], t['ctotal_g'], 10)
            ct += 1
        self.assertTrue(ct>0)

        if os.path.exists(baselineCatName):
            os.unlink(baselineCatName)

        if os.path.exists(testCatName):
            os.unlink(testCatName)

    def testPhotometricIndicesRaw(self):
        """
        Use manMagCalc_list with specified indices on an Sed.  Make sure
        that the appropriate magnitudes are or are not Nan
        """
        starName = os.path.join(lsst.utils.getPackageDir('sims_sed_library'),defaultSpecMap['km20_5750.fits_g40_5790'])
        starPhot = loadTotalBandpassesFromFiles()
        testSed = Sed()
        testSed.readSED_flambda(starName)
        indices = [1,3]
        mags = starPhot.calcMagList(testSed, indices=indices)
        self.assertTrue(numpy.isnan(mags[0]))
        self.assertFalse(numpy.isnan(mags[1]))
        self.assertTrue(numpy.isnan(mags[2]))
        self.assertFalse(numpy.isnan(mags[3]))
        self.assertTrue(numpy.isnan(mags[4]))
        self.assertTrue(numpy.isnan(mags[5]))
        self.assertTrue(len(mags)==6)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(variabilityUnitTest)
    suites += unittest.makeSuite(photometryUnitTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
