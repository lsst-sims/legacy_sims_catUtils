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
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import PhotometricParameters, calcSNR_m5, LSSTdefaults
from lsst.sims.photUtils import calcSNR_sed, magErrorFromSNR
from lsst.sims.photUtils.utils import setM5
from lsst.sims.catUtils.utils import cartoonStars, cartoonGalaxies, testStars, testGalaxies, \
                                     cartoonStarsOnlyI, cartoonStarsIZ, \
                                     cartoonGalaxiesIG, galaxiesWithHoles
from lsst.sims.catUtils.mixins import PhotometryBase, PhotometryStars, PhotometryGalaxies


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

        cartoonDir = os.path.join(lsst.utils.getPackageDir('sims_photUtils'), 'tests', 'cartoonSedTestData')
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
        self.assertTrue(len(mags)==len(test_cat.cartoonBandpassDict))
        self.assertTrue(len(mags)>0)
        for j in range(len(mags)):
            self.assertAlmostEqual(mags[j],test_cat.magnitudeMasterList[i][j],10)
        i += 1

        os.unlink("testStarsCartoon.txt")

    def testAlternateBandpassesGalaxies(self):
        """
        the same as testAlternateBandpassesStars, but for galaxies
        """

        catName = 'testAlternateBandpassesGalaxies.txt'

        obs_metadata_pointed=ObservationMetaData(mjd=50000.0,
                               boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                               boundLength=10.0)

        test_cat=cartoonGalaxies(self.galaxy,obs_metadata=obs_metadata_pointed)
        test_cat.write_catalog(catName)

        dtype = numpy.dtype([
                             ('galid', numpy.int),
                             ('ra', numpy.float),
                             ('dec', numpy.float),
                             ('uTotal', numpy.float),
                             ('gTotal', numpy.float),
                             ('rTotal', numpy.float),
                             ('iTotal', numpy.float),
                             ('zTotal', numpy.float),
                             ('uBulge', numpy.float),
                             ('gBulge', numpy.float),
                             ('rBulge', numpy.float),
                             ('iBulge', numpy.float),
                             ('zBulge', numpy.float),
                             ('uDisk', numpy.float),
                             ('gDisk', numpy.float),
                             ('rDisk', numpy.float),
                             ('iDisk', numpy.float),
                             ('zDisk', numpy.float),
                             ('uAgn', numpy.float),
                             ('gAgn', numpy.float),
                             ('rAgn', numpy.float),
                             ('iAgn', numpy.float),
                             ('zAgn', numpy.float),
                             ('bulgeName', str, 200),
                             ('bulgeNorm', numpy.float),
                             ('bulgeAv', numpy.float),
                             ('diskName', str, 200),
                             ('diskNorm', numpy.float),
                             ('diskAv', numpy.float),
                             ('agnName', str, 200),
                             ('agnNorm', numpy.float),
                             ('redshift', numpy.float)
                            ])

        catData = numpy.genfromtxt(catName, dtype=dtype, delimiter=', ')

        cartoonDir = lsst.utils.getPackageDir('sims_photUtils')
        cartoonDir = os.path.join(cartoonDir, 'tests', 'cartoonSedTestData')
        sedDir = lsst.utils.getPackageDir('sims_sed_library')

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
                    self.assertTrue(numpy.isnan(line['%sBulge' % bp]))
                    bulgeMagList.append(numpy.NaN)
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
                    self.assertTrue(numpy.isnan(line['%sDisk' % bp]))
                    diskMagList.append(numpy.NaN)
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
                    self.assertTrue(numpy.isnan(line['%sAgn' % bp]))
                    agnMagList.append(numpy.NaN)
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

            totalMags = PhotometryGalaxies().sum_magnitudes(bulge=numpy.array(bulgeMagList),
                                                            disk=numpy.array(diskMagList),
                                                            agn=numpy.array(agnMagList))

            for testMag, bpName in zip(totalMags, keys):
                if numpy.isnan(line['%sTotal' % bpName]):
                    self.assertTrue(numpy.isnan(testMag))
                else:
                    self.assertAlmostEqual(testMag, line['%sTotal' % bpName],10)


        self.assertTrue(ct>0)
        if os.path.exists(catName):
            os.unlink(catName)

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
        starPhot = BandpassDict.loadTotalBandpassesFromFiles()
        testSed = Sed()
        testSed.readSED_flambda(starName)
        indices = [1,3]
        mags = starPhot.magListForSed(testSed, indices=indices)
        self.assertTrue(numpy.isnan(mags[0]))
        self.assertFalse(numpy.isnan(mags[1]))
        self.assertTrue(numpy.isnan(mags[2]))
        self.assertFalse(numpy.isnan(mags[3]))
        self.assertTrue(numpy.isnan(mags[4]))
        self.assertTrue(numpy.isnan(mags[5]))
        self.assertTrue(len(mags)==6)


class UncertaintyMixinTest(unittest.TestCase):
    def setUp(self):
        starName = os.path.join(lsst.utils.getPackageDir('sims_sed_library'),defaultSpecMap['km20_5750.fits_g40_5790'])
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
            hardwareComponents.append(os.path.join(lsst.utils.getPackageDir('throughputs'),'baseline',c))

        self.bandpasses = ['u', 'g', 'r', 'i', 'z', 'y']
        for b in self.bandpasses:
            filterName = os.path.join(lsst.utils.getPackageDir('throughputs'),'baseline','filter_%s.dat' % b)
            components = hardwareComponents + [filterName]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.hardwareBandpasses.append(bandpassDummy)
            components = components + [os.path.join(lsst.utils.getPackageDir('throughputs'),'baseline','atmos.dat')]
            bandpassDummy = Bandpass()
            bandpassDummy.readThroughputList(components)
            self.totalBandpasses.append(bandpassDummy)


    def testUncertaintyExceptions(self):
        """
        Test the calculateMagnitudeUncertainty raises exceptions when it needs to
        """
        phot = PhotometryBase()
        totalDict, hardwareDict = BandpassDict.loadBandpassesFromFiles()
        magnitudes = numpy.array([22.0, 23.0, 24.0, 25.0, 26.0, 27.0])
        shortMagnitudes = numpy.array([22.0])
        self.assertRaises(RuntimeError, phot.calculateMagnitudeUncertainty, magnitudes, totalDict)
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, bandpassName='g', m5=23.0)
        self.assertRaises(RuntimeError, phot.calculateMagnitudeUncertainty, shortMagnitudes, totalDict, \
                          obs_metadata=obs_metadata)

        photParams = PhotometricParameters()
        shortGamma = numpy.array([1.0, 1.0])
        self.assertRaises(RuntimeError, calcSNR_m5, magnitudes, totalDict.values(), shortMagnitudes, photParams)
        self.assertRaises(RuntimeError, calcSNR_m5, shortMagnitudes, totalDict.values(), magnitudes, photParams)
        self.assertRaises(RuntimeError, calcSNR_m5, magnitudes, totalDict.values(), magnitudes, photParams, gamma=shortGamma)
        snr, gg = calcSNR_m5(magnitudes, totalDict.values(), magnitudes, photParams)




    def testRawUncertainty(self):
        """
        Test that values calculated by calculatePhotometricUncertainty agree
        with values calculated by calcSNR_sed
        """

        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        bandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
        phot = PhotometryBase()
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.bandpasses)
        magnitudes = bandpassDict.magListForSed(self.starSED)

        skySeds = []

        for i in range(len(self.bandpasses)):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[self.bandpasses[i]], skyDummy,
                                       self.totalBandpasses[i], self.hardwareBandpasses[i],
                                       seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                                       photParams=PhotometricParameters())

            skySeds.append(normalizedSkyDummy)

        sigma = phot.calculateMagnitudeUncertainty(magnitudes, bandpassDict, obs_metadata=obs_metadata)
        for i in range(len(self.bandpasses)):
            snr = calcSNR_sed(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            ss = 2.5*numpy.log10(1.0+1.0/snr)
            ss = numpy.sqrt(ss*ss + numpy.power(phot.photParams.sigmaSys,2))
            msg = '%e is not %e; failed' % (ss, sigma[i])
            self.assertAlmostEqual(ss, sigma[i], 10, msg=msg)

    def testSystematicUncertainty(self):
        """
        Test that systematic uncertainty is added correctly.
        """
        sigmaSys = 0.002
        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        photParams= PhotometricParameters(sigmaSys=sigmaSys)

        phot = PhotometryBase()
        phot.photParams = photParams

        bandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.bandpasses)
        magnitudes = bandpassDict.magListForSed(self.starSED)

        skySeds = []

        for i in range(len(self.bandpasses)):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[self.bandpasses[i]], skyDummy,
                                                       self.totalBandpasses[i], self.hardwareBandpasses[i],
                                                       seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                                                       photParams=PhotometricParameters())

            skySeds.append(normalizedSkyDummy)

        sigma = phot.calculateMagnitudeUncertainty(magnitudes, bandpassDict, obs_metadata=obs_metadata)
        for i in range(len(self.bandpasses)):
            snr = calcSNR_sed(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            testSNR, gamma = calcSNR_m5(numpy.array([magnitudes[i]]), [self.totalBandpasses[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters(sigmaSys=0.0))

            self.assertAlmostEqual(snr, testSNR[0], 10, msg = 'failed on calcSNR_m5 test %e != %e ' \
                                                               % (snr, testSNR[0]))

            control = numpy.sqrt(numpy.power(magErrorFromSNR(testSNR),2) + numpy.power(sigmaSys,2))

            msg = '%e is not %e; failed' % (sigma[i], control)

            self.assertAlmostEqual(sigma[i], control, 10, msg=msg)


    def testNoSystematicUncertainty(self):
        """
        Test that systematic uncertainty is handled correctly when set to None.
        """
        m5 = [23.5, 24.3, 22.1, 20.0, 19.5, 21.7]
        photParams= PhotometricParameters(sigmaSys=0.0)

        phot = PhotometryBase()
        phot.photParams = photParams

        bandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
        obs_metadata = ObservationMetaData(unrefractedRA=23.0, unrefractedDec=45.0, m5=m5, bandpassName=self.bandpasses)
        magnitudes = bandpassDict.magListForSed(self.starSED)

        skySeds = []

        for i in range(len(self.bandpasses)):
            skyDummy = Sed()
            skyDummy.readSED_flambda(os.path.join(lsst.utils.getPackageDir('throughputs'), 'baseline', 'darksky.dat'))
            normalizedSkyDummy = setM5(obs_metadata.m5[self.bandpasses[i]], skyDummy,
                                                       self.totalBandpasses[i], self.hardwareBandpasses[i],
                                                       seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                                                       photParams=PhotometricParameters())

            skySeds.append(normalizedSkyDummy)

        sigma = phot.calculateMagnitudeUncertainty(magnitudes, bandpassDict, obs_metadata=obs_metadata)
        for i in range(len(self.bandpasses)):
            snr = calcSNR_sed(self.starSED, self.totalBandpasses[i], skySeds[i], self.hardwareBandpasses[i],
                              seeing=LSSTdefaults().seeing(self.bandpasses[i]),
                              photParams=PhotometricParameters())

            testSNR, gamma = calcSNR_m5(numpy.array([magnitudes[i]]), [self.totalBandpasses[i]],
                                           numpy.array([m5[i]]), photParams=PhotometricParameters(sigmaSys=0.0))

            self.assertAlmostEqual(snr, testSNR[0], 10, msg = 'failed on calcSNR_m5 test %e != %e ' \
                                                               % (snr, testSNR[0]))

            control = magErrorFromSNR(testSNR)

            msg = '%e is not %e; failed' % (sigma[i], control)

            self.assertAlmostEqual(sigma[i], control, 10, msg=msg)



def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(variabilityUnitTest)
    suites += unittest.makeSuite(photometryUnitTest)
    suites += unittest.makeSuite(UncertaintyMixinTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
