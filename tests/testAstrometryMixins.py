import numpy

import os
import unittest
import warnings
import sys
import math
import palpy as pal
import lsst.utils.tests as utilsTests

from lsst.utils import getPackageDir
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.utils import ObservationMetaData, arcsecFromRadians
from lsst.sims.utils import _observedFromAppGeo, _pupilCoordsFromRaDec
from lsst.sims.coordUtils import focalPlaneCoordsFromPupilCoords, _focalPlaneCoordsFromRaDec
from lsst.sims.coordUtils import chipNameFromPupilCoords, _chipNameFromRaDec
from lsst.sims.coordUtils import pixelCoordsFromPupilCoords, _pixelCoordsFromRaDec
from lsst.sims.catalogs.generation.utils import myTestStars, makeStarTestDB, \
                                                myTestGals, makeGalTestDB
import lsst.afw.cameraGeom.testUtils as camTestUtils
from lsst.sims.catUtils.mixins import AstrometryStars, AstrometryGalaxies, CameraCoords


class AstrometryTestStars(myTestStars):
    database = 'AstrometryTestStarDatabase.db'

class AstrometryTestGalaxies(myTestGals):
    database = 'AstrometryTestGalaxyDatabase.db'

class parallaxTestCatalog(InstanceCatalog, AstrometryStars):
    column_outputs = ['raJ2000', 'decJ2000', 'raObserved', 'decObserved',
                      'properMotionRa', 'properMotionDec',
                      'radialVelocity', 'parallax']

    transformations = {'raJ2000':numpy.degrees, 'decJ2000':numpy.degrees,
                       'raObserved':numpy.degrees, 'decObserved':numpy.degrees,
                       'properMotionRa':numpy.degrees, 'properMotionDec':numpy.degrees,
                       'parallax':arcsecFromRadians}

    default_formats = {'f':'%.12f'}

class testCatalog(InstanceCatalog,AstrometryStars,CameraCoords):
    """
    A (somewhat meaningless) instance catalog class that will allow us
    to run the astrometry routines for testing purposes
    """
    catalog_type = 'test_stars'
    column_outputs=['id','raPhoSim','decPhoSim','raObserved','decObserved',
                   'x_pupil','y_pupil',
                   'chipName', 'xPix', 'yPix','xFocalPlane','yFocalPlane']
    #Needed to do camera coordinate transforms.
    camera = camTestUtils.CameraWrapper().camera
    default_formats = {'f':'%.12f'}

    delimiter = ';' #so that numpy.loadtxt can parse the chipNames which may contain commas
                     #(see testClassMethods)

    default_columns = [('properMotionRa', 0., float),
                       ('properMotionDec', 0., float),
                       ('parallax', 1.2, float),
                       ('radial_velocity', 0., float)]


class testStellarCatalog(InstanceCatalog, AstrometryStars, CameraCoords):
    """
    Define a catalog of stars with all possible astrometric columns
    """

    camera = camTestUtils.CameraWrapper().camera

    column_outputs = ['glon', 'glat',
                      'x_pupil', 'y_pupil',
                      'xPix', 'yPix',
                      'xFocalPlane', 'yFocalPlane',
                      'chipName', 'raPhoSim', 'decPhoSim',
                      'raObserved', 'decObserved']

class testGalaxyCatalog(InstanceCatalog, AstrometryGalaxies, CameraCoords):
    """
    Define a catalog of galaxies with all possible astrometric columns
    """

    camera = camTestUtils.CameraWrapper().camera

    column_outputs = ['glon', 'glat',
                      'x_pupil', 'y_pupil', 'xPix', 'yPix', 'xFocalPlane', 'yFocalPlane',
                      'chipName', 'raPhoSim', 'decPhoSim', 'raObserved', 'decObserved']

class astrometryUnitTest(unittest.TestCase):
    """
    The bulk of this unit test involves inputting a set list of input values
    and comparing the astrometric results to results derived from SLALIB run
    with the same input values.  We have to create a test catalog artificially (rather than
    querying the database) because SLALIB was originally run on values that did not correspond
    to any particular Opsim run.
    """

    @classmethod
    def setUpClass(cls):
        # Create test databases
        cls.starDBName = 'AstrometryTestStarDatabase.db'
        cls.galDBName = 'AstrometryTestGalaxyDatabase.db'
        if os.path.exists(cls.starDBName):
            os.unlink(cls.starDBName)
        makeStarTestDB(filename=cls.starDBName,
                      size=100000, seedVal=1, ramin=199.98*math.pi/180., dra=0.04*math.pi/180.)

        if os.path.exists(cls.galDBName):
            os.unlink(cls.galDBName)
        makeGalTestDB(filename=cls.galDBName,
                      size=100000, seedVal=1, ramin=199.98*math.pi/180., dra=0.04*math.pi/180.)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.starDBName):
            os.unlink(cls.starDBName)

        if os.path.exists(cls.galDBName):
            os.unlink(cls.galDBName)

    def setUp(self):
        self.starDBObject = AstrometryTestStars()
        self.galaxyDBObject = AstrometryTestGalaxies()
        self.metadata={}

        #below are metadata values that need to be set in order for
        #get_getFocalPlaneCoordinates to work.  If we had been querying the database,
        #these would be set to meaningful values.  Because we are generating
        #an artificial set of inputs that must comport to the baseline SLALIB
        #inputs, these are set arbitrarily by hand
        self.metadata['pointingRA'] = (numpy.radians(200.0), float)
        self.metadata['pointingDec'] = (numpy.radians(-30.0), float)
        self.metadata['Opsim_rotskypos'] = (1.0, float)

        self.obs_metadata=ObservationMetaData(mjd=57388.0,
                                     boundType='circle',
                                     boundLength=0.05,
                                     phoSimMetaData=self.metadata)

        self.cat = testCatalog(self.starDBObject, obs_metadata=self.obs_metadata)
        self.tol=1.0e-5

    @classmethod
    def tearDownClass(cls):
        if os.path.exists('AstrometryTestDatabase.db'):
            os.unlink('AstrometryTestDatabase.db')

    def tearDown(self):
        del self.starDBObject
        del self.galaxyDBObject
        del self.cat
        del self.obs_metadata
        del self.metadata
        del self.tol


    def testWritingOfStars(self):
        """
        Try writing a catalog with all possible Astrometric columns
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'starWritingTest.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        stars = testStellarCatalog(self.starDBObject, obs_metadata=self.obs_metadata)

        stars.write_catalog(catName)

        dtypeList = [(name, numpy.float) for name in stars._column_outputs]
        testData = numpy.genfromtxt(catName, delimiter=', ', dtype=numpy.dtype(dtypeList))

        self.assertGreater(len(testData), 0)

        if os.path.exists(catName):
            os.unlink(catName)

    def testWritingOfGalaxies(self):
        """
        Try writing a catalog with all possible Astrometric columns
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'galWritingTest.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        galaxies = testGalaxyCatalog(self.galaxyDBObject, obs_metadata=self.obs_metadata)
        galaxies.write_catalog(catName)

        dtypeList = [(name, numpy.float) for name in galaxies._column_outputs]
        testData = numpy.genfromtxt(catName, dtype=numpy.dtype(dtypeList))

        self.assertGreater(len(testData), 0)

        if os.path.exists(catName):
            os.unlink(catName)



    def testUtilityMethods(self):
        """
        Generate a catalog using the methods from AstrometryUtils.py and CameraUtils.py.
        Read that data in, and then recalculate the values 'by hand' to make sure
        that they are consistent.
        """

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'AstrometryUtilityCatalog.txt')

        if os.path.exists(catName):
            os.unlink(catName)

        self.cat.write_catalog(catName)

        dtype = [('id',int),
                 ('raPhoSim',float), ('decPhoSim',float),
                 ('raObserved',float), ('decObserved',float),
                 ('x_pupil',float), ('y_pupil',float), ('chipName',str,11),
                 ('xPix',float), ('yPix',float),
                 ('xFocalPlane',float), ('yFocalPlane',float)]

        baselineData = numpy.loadtxt(catName, dtype=dtype, delimiter=';')

        self.assertGreater(len(baselineData), 0)

        pupilTest = _pupilCoordsFromRaDec(baselineData['raObserved'],
                                              baselineData['decObserved'],
                                              obs_metadata=self.obs_metadata,
                                              epoch=2000.0)

        for (xxtest, yytest, xx, yy) in \
                zip(pupilTest[0], pupilTest[1], baselineData['x_pupil'], baselineData['y_pupil']):
            self.assertAlmostEqual(xxtest,xx,6)
            self.assertAlmostEqual(yytest,yy,6)

        focalTest = focalPlaneCoordsFromPupilCoords(pupilTest[0], pupilTest[1], camera=self.cat.camera)

        focalRa = _focalPlaneCoordsFromRaDec(baselineData['raObserved'], baselineData['decObserved'],
                                             epoch=self.cat.db_obj.epoch, obs_metadata=self.cat.obs_metadata,
                                             camera=self.cat.camera)

        for (xxtest, yytest, xxra, yyra, xx, yy) in \
                zip(focalTest[0], focalTest[1], focalRa[0], focalRa[1],
                        baselineData['xFocalPlane'], baselineData['yFocalPlane']):

            self.assertAlmostEqual(xxtest,xx,6)
            self.assertAlmostEqual(yytest,yy,6)
            self.assertAlmostEqual(xxra,xx,6)
            self.assertAlmostEqual(yyra,yy,6)

        pixTest = pixelCoordsFromPupilCoords(pupilTest[0], pupilTest[1], camera=self.cat.camera)

        pixTestRaDec = _pixelCoordsFromRaDec(baselineData['raObserved'], baselineData['decObserved'],
                                             epoch=self.cat.db_obj.epoch,
                                             obs_metadata=self.cat.obs_metadata,
                                             camera=self.cat.camera)

        for (xxtest, yytest, xxra, yyra, xx, yy) in \
                zip(pixTest[0], pixTest[1], pixTestRaDec[0], pixTestRaDec[1],
                           baselineData['xPix'], baselineData['yPix']):

            if not numpy.isnan(xx) and not numpy.isnan(yy):
                self.assertAlmostEqual(xxtest,xx,5)
                self.assertAlmostEqual(yytest,yy,5)
                self.assertAlmostEqual(xxra,xx,5)
                self.assertAlmostEqual(yyra,yy,5)
            else:
                self.assertTrue(numpy.isnan(xx))
                self.assertTrue(numpy.isnan(yy))
                self.assertTrue(numpy.isnan(xxra))
                self.assertTrue(numpy.isnan(yyra))
                self.assertTrue(numpy.isnan(xxtest))
                self.assertTrue(numpy.isnan(yytest))


        nameTest = chipNameFromPupilCoords(pupilTest[0], pupilTest[1],
                                           camera=self.cat.camera)

        nameRA = _chipNameFromRaDec(baselineData['raObserved'], baselineData['decObserved'],
                                    epoch=self.cat.db_obj.epoch, obs_metadata=self.cat.obs_metadata,
                                    camera=self.cat.camera)

        for (ntest, nra, ncontrol) in zip(nameTest, nameRA, baselineData['chipName']):
            if ncontrol != 'None':
                self.assertEqual(ntest,ncontrol)
                self.assertEqual(nra,ncontrol)
            else:
                self.assertTrue(ntest is None)
                self.assertTrue(nra is None)

        if os.path.exists(catName):
            os.unlink(catName)


    def testParallax(self):
        """
        This test will output a catalog of ICRS and observed positions.
        It will also output the quantities (proper motion, radial velocity,
        and parallax) needed to apply the transformaiton between the two.
        It will then run the catalog through PALPY and verify that the catalog
        generating code correctly applied the transformations.
        """

        #create and write a catalog that performs astrometric transformations
        #on a cartoon star database
        cat = parallaxTestCatalog(self.starDBObject, obs_metadata=self.obs_metadata)
        parallaxName = os.path.join(getPackageDir('sims_catUtils'), 'tests',
                                    'scratchSpace', 'parallaxCatalog.sav')

        if os.path.exists(parallaxName):
            os.unlink(parallaxName)

        cat.write_catalog(parallaxName)

        data = numpy.genfromtxt(parallaxName,delimiter=',')
        self.assertGreater(len(data), 0)

        epoch = cat.db_obj.epoch
        mjd = cat.obs_metadata.mjd
        prms = pal.mappa(epoch, mjd)
        for vv in data:
            #run the PALPY routines that actuall do astrometry `by hand' and compare
            #the results to the contents of the catalog
            ra0 = numpy.radians(vv[0])
            dec0 = numpy.radians(vv[1])
            pmra = numpy.radians(vv[4])
            pmdec = numpy.radians(vv[5])
            rv = vv[6]
            px = vv[7]
            ra_apparent, dec_apparent = pal.mapqk(ra0, dec0, pmra, pmdec, px, rv, prms)
            ra_apparent = numpy.array([ra_apparent])
            dec_apparent = numpy.array([dec_apparent])
            raObserved, decObserved = _observedFromAppGeo(ra_apparent, dec_apparent,
                                                                 obs_metadata=cat.obs_metadata)

            self.assertAlmostEqual(raObserved[0],numpy.radians(vv[2]),7)
            self.assertAlmostEqual(decObserved[0],numpy.radians(vv[3]),7)

        if os.path.exists(parallaxName):
            os.unlink(parallaxName)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(astrometryUnitTest)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
