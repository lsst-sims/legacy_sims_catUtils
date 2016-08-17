import numpy as np

import os
import unittest
import lsst
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import AstrometryStars, AstrometryGalaxies
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies
from lsst.sims.catUtils.utils import setupPhotometryCatalog
from lsst.sims.catUtils.utils import makeStarDatabase, makeGalaxyDatabase


def setup_module(module):
    lsst.utils.tests.init()


class testStarCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):
    """
    A class with no photometry columns.  Meant to be passed to setupPhotometryCatalog
    where it will be given photometry columns
    """
    column_outputs = ['raObserved', 'decObserved']
    default_formats = {'f':'%.12e'}

class baselineStarCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):
    """
    Baseline photometry catalog against which to compare testStarCatalog
    """
    column_outputs = ['raObserved', 'decObserved']
    default_formats = {'f':'%.12e'}

class testGalaxyCatalog(InstanceCatalog, AstrometryGalaxies, PhotometryGalaxies):
    """
    A class with no photometry columns.  Meant to be passed to setupPhotometryCatalog
    where it will be given photometry columns
    """
    column_outputs = ['raObserved', 'decObserved']
    default_formats = {'f':'%.12e'}

class baselineGalaxyCatalog(InstanceCatalog, AstrometryGalaxies, PhotometryGalaxies):
    """
    Baseline photometry catalog against which to compare testGalaxyCatalog
    """
    column_outputs = ['raObserved', 'decObserved']
    default_formats = {'f':'%.12e'}

class testStarDBObject(CatalogDBObject):
    """
    CatalogDBObject to map our test database of stars
    """
    tableid = 'StarAllForceseek'
    idColKey = 'id'
    raColName = 'ra'
    decColName = 'decl'
    objectTypeId = 49
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('magNorm', None),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', '3.1*ebv'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40)]

class testGalaxyDBObject(CatalogDBObject):

    #: This is the base table for the galaxies
    #tableid = 'final_clone_db'
    tableid = 'galaxy'
    idColKey = 'galtileid'
    raColName = '((CAST(ra AS NUMERIC(9,6))%360.)+360.)%360.'
    decColName = 'dec'
    objectTypeId = 51

    columns = [('galtileid', None, np.int64),
            ('galid', None, str, 30),
            ('raJ2000', 'ra*PI()/180.'),
            ('decJ2000', 'dec*PI()/180.'),
            ('raJ2000Bulge', 'bra*PI()/180.'),
            ('decJ2000Bulge', 'bdec*PI()/180.'),
            ('raJ2000Disk', 'dra*PI()/180.'),
            ('decJ2000Disk', 'ddec*PI()/180.'),
            ('raJ2000Agn', 'agnra*PI()/180.'),
            ('decJ2000Agn', 'agndec*PI()/180.'),
            ('magNormBulge', 'magnorm_bulge'),
            ('magNormDisk', 'magnorm_disk'),
            ('magNormAgn', 'magnorm_agn'),
            ('sedFilenameBulge', 'sedname_bulge', unicode, 40),
            ('sedFilenameDisk', 'sedname_disk', unicode, 40),
            ('sedFilenameAgn', 'sedname_agn', unicode, 40),
            ('majorAxisBulge', 'a_b*PI()/648000.'),
            ('minorAxisBulge', 'b_b*PI()/648000.'),
            ('positionAngleBulge', 'pa_bulge*PI()/180.'),
            ('sindexBulge', 'bulge_n', int),
            ('majorAxisDisk', 'a_d*PI()/648000.'),
            ('minorAxisDisk', 'b_d*PI()/648000.'),
            ('positionAngleDisk', 'pa_disk*PI()/180.'),
            ('sindexDisk', 'disk_n', int),
            ('internalExtinctionModelBulge', 'ext_model_b', str, 3),
            ('internalAvBulge', 'av_b'),
            ('internalRvBulge', 'rv_b'),
            ('internalExtinctionModelDisk', 'ext_model_d', str, 3),
            ('internalAvDisk', 'av_d'),
            ('internalRvDisk', 'rv_d'),
            ('lsst_u', 'u_ab'),
            ('lsst_g', 'g_ab'),
            ('lsst_r', 'r_ab'),
            ('lsst_i', 'i_ab'),
            ('lsst_z', 'z_ab'),
            ('lsst_y', 'y_ab')]

class InstanceCatalogSetupUnittest(unittest.TestCase):

    def setUp(self):
        self.driver = 'sqlite'
        self.StarDBName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                      'testSetup_setupTestStars.db')

        self.GalaxyDBName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                         'testSetup_setupTestGalaxies.db')

        if os.path.exists(self.StarDBName):
            os.unlink(self.StarDBName)

        if os.path.exists(self.GalaxyDBName):
            os.unlink(self.GalaxyDBName)

        self.pointingRA = 50.0
        self.pointingDec = -5.0
        self.radius = 1.0
        makeStarDatabase(filename=self.StarDBName, size=100,
                         pointingRA=self.pointingRA,
                         pointingDec=self.pointingDec,
                         radius=self.radius)


        makeGalaxyDatabase(filename=self.GalaxyDBName, size=100,
                           pointingRA=self.pointingRA,
                           pointingDec=self.pointingDec,
                           radius=self.radius)

        self.starDBObj = testStarDBObject(driver=self.driver, database= self.StarDBName)
        self.galaxyDBObj = testGalaxyDBObject(driver=self.driver, database=self.GalaxyDBName)

        self.obs_metadata = ObservationMetaData(pointingRA=self.pointingRA,
                                                pointingDec=self.pointingDec,
                                                boundType='circle', boundLength=self.radius,
                                                bandpassName='g', mjd=57000.0,
                                                m5=24.5)

        self.obs_metadata_compound = ObservationMetaData(pointingRA=self.pointingRA,
                                                         pointingDec=self.pointingDec,
                                                         boundType='circle', boundLength=self.radius,
                                                         bandpassName=['g','i'], mjd=57000.0,
                                                         m5=[24.5, 17.5])

    def tearDown(self):
        if os.path.exists(self.StarDBName):
            os.unlink(self.StarDBName)

        if os.path.exists(self.GalaxyDBName):
            os.unlink(self.GalaxyDBName)

        del self.starDBObj
        del self.galaxyDBObj
        del self.StarDBName
        del self.GalaxyDBName
        del self.pointingRA
        del self.pointingDec
        del self.radius
        del self.obs_metadata

    def testExceptions(self):
        """
        Make sure that setupPhotometryCatalog throws errors when it is supposed to
        """

        class dummyClass(object):
            def __init__(self):
                pass

        xx = dummyClass()
        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=xx,
                          dbConnection=self.starDBObj, catalogClass=testStarCatalog)

        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=self.obs_metadata,
                          dbConnection=xx, catalogClass=testStarCatalog)

        self.assertRaises(RuntimeError, setupPhotometryCatalog, obs_metadata=self.obs_metadata,
                          dbConnection=self.starDBObj, catalogClass=dummyClass)


    def testSetupPhotometry(self):
        """
        Make sure that catalogs instantiated by setupPhotometryCatalog contain the
        correct columns.
        """

        #test case with a single bandpass
        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata, dbConnection=self.starDBObj,
                                     catalogClass=testStarCatalog)

        self.assertIn('lsst_g', cat.iter_column_names())
        self.assertNotIn('lsst_u', cat.iter_column_names())
        self.assertNotIn('lsst_r', cat.iter_column_names())
        self.assertNotIn('lsst_i', cat.iter_column_names())
        self.assertNotIn('lsst_z', cat.iter_column_names())
        self.assertNotIn('lsst_y', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_g', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_u', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_r', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_i', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_z', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_y', cat.iter_column_names())

        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata, dbConnection=self.starDBObj,
                                     catalogClass=testStarCatalog, uncertainty=True)

        self.assertIn('lsst_g', cat.iter_column_names())
        self.assertNotIn('lsst_u', cat.iter_column_names())
        self.assertNotIn('lsst_r', cat.iter_column_names())
        self.assertNotIn('lsst_i', cat.iter_column_names())
        self.assertNotIn('lsst_z', cat.iter_column_names())
        self.assertNotIn('lsst_y', cat.iter_column_names())
        self.assertIn('sigma_lsst_g', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_u', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_r', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_i', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_z', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_y', cat.iter_column_names())

        #test case with two bandpasses
        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                     dbConnection=self.starDBObj, catalogClass=testStarCatalog)

        self.assertIn('lsst_g', cat.iter_column_names())
        self.assertIn('lsst_i', cat.iter_column_names())
        self.assertNotIn('lsst_u', cat.iter_column_names())
        self.assertNotIn('lsst_r', cat.iter_column_names())
        self.assertNotIn('lsst_z', cat.iter_column_names())
        self.assertNotIn('lsst_y', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_g', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_u', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_r', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_i', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_z', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_y', cat.iter_column_names())

        cat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                     dbConnection=self.starDBObj, catalogClass=testStarCatalog,
                                     uncertainty=True)

        self.assertIn('lsst_g', cat.iter_column_names())
        self.assertIn('lsst_i', cat.iter_column_names())
        self.assertNotIn('lsst_u', cat.iter_column_names())
        self.assertNotIn('lsst_r', cat.iter_column_names())
        self.assertNotIn('lsst_z', cat.iter_column_names())
        self.assertNotIn('lsst_y', cat.iter_column_names())
        self.assertIn('sigma_lsst_g', cat.iter_column_names())
        self.assertIn('sigma_lsst_i', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_u', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_r', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_z', cat.iter_column_names())
        self.assertNotIn('sigma_lsst_y', cat.iter_column_names())

        #make sure that class default columns did not get overwritten
        cat = testStarCatalog(self.starDBObj, obs_metadata=self.obs_metadata)

        self.assertNotIn('lsst_u', cat.iter_column_names())
        self.assertNotIn('lsst_g', cat.iter_column_names())
        self.assertNotIn('lsst_r', cat.iter_column_names())
        self.assertNotIn('lsst_i', cat.iter_column_names())
        self.assertNotIn('lsst_z', cat.iter_column_names())
        self.assertNotIn('lsst_y', cat.iter_column_names())



    def testActualCatalog(self):
        """
        Make sure that the values written to catalogs that are instantiated using
        setupPhotometryCatalog are correct
        """
        msgroot = ['failed on stars; ', 'failed on galaxies; ']

        testCatClasses = [testStarCatalog, testGalaxyCatalog]
        testCatDBs = [self.starDBObj, self.galaxyDBObj]
        baselineCats = []
        baselineCats.append(baselineStarCatalog(self.starDBObj, obs_metadata=self.obs_metadata_compound,
                                                column_outputs=['lsst_g', 'sigma_lsst_g', 'lsst_i', 'sigma_lsst_i']))

        baselineCats.append(baselineGalaxyCatalog(self.galaxyDBObj, obs_metadata=self.obs_metadata_compound,
                                                  column_outputs=['lsst_g', 'sigma_lsst_g', 'lsst_i', 'sigma_lsst_i']))

        testName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'testSetUp_testActual_testSetupCat.txt')

        baseName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                'testSetUp_testActual_baseSetupCat.txt')

        if os.path.exists(testName):
            os.unlink(testName)

        if os.path.exists(baseName):
            os.unlink(baseName)


        basedtype = np.dtype([('raObserved', np.float), ('decObserved', np.float),
                                 ('lsst_g', np.float), ('sigma_lsst_g', np.float),
                                 ('lsst_i', np.float), ('sigma_lsst_i', np.float)])

        for (testCatClass, dbo, baselineCat, msgr) in zip(testCatClasses, testCatDBs, baselineCats, msgroot):

            testdtype = np.dtype([('raObserved', np.float), ('decObserved', np.float),
                                 ('lsst_g', np.float)])


            testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata,
                                              dbConnection=dbo,
                                              catalogClass=testCatClass)

            testCat.write_catalog(testName)
            baselineCat.write_catalog(baseName)

            testData = np.genfromtxt(testName, dtype=testdtype, delimiter=',')
            baseData = np.genfromtxt(baseName, dtype=basedtype, delimiter=',')
            self.assertGreater(len(testData), 0)
            self.assertGreater(len(baseData), 0)

            ct = 0
            for b, t in zip(baseData, testData):
                self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12,
                                       msg = '%s single column; %.12e != %.12e' % (msgr, b['lsst_g'], t['lsst_g']))
                ct +=1

            self.assertGreater(ct, 0)

            testdtype = np.dtype([('raObserved', np.float), ('decObserved', np.float),
                                     ('lsst_g', np.float), ('lsst_i', np.float)])

            testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                             dbConnection=dbo,
                                             catalogClass=testCatClass)
            testCat.write_catalog(testName)
            testData = np.genfromtxt(testName, dtype=testdtype, delimiter=',')
            self.assertGreater(len(testData), 0)
            ct = 0
            for b, t in zip(baseData, testData):
                self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['lsst_g'], t['lsst_g']))
                self.assertAlmostEqual(b['lsst_i'], t['lsst_i'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['lsst_i'], t['lsst_i']))
                ct += 1

            self.assertGreater(ct, 0)

            if os.path.exists(testName):
                os.unlink(testName)
            if os.path.exists(baseName):
                os.unlink(baseName)

    def testActualCatalogWithUncertainty(self):
        """
        Make sure that the values written to catalogs that are instantiated using
        setupPhotometryCatalog are correct (include photometric uncertainty)
        """

        msgroot = ['failed on stars; ', 'failed on galaxies; ']

        testCatClasses = [testStarCatalog, testGalaxyCatalog]
        testCatDBs = [self.starDBObj, self.galaxyDBObj]
        baselineCats = []

        #need to set up the baseline catalogs with the compound obs_metadata so that they get the
        #correct m5 values for both magnitudes (otherwise, they will use LSST defaults, which
        #disagree with our cartoon test case)
        baselineCats.append(baselineStarCatalog(self.starDBObj, obs_metadata=self.obs_metadata_compound,
                                                column_outputs=['lsst_g', 'lsst_i', 'sigma_lsst_g', 'sigma_lsst_i']))

        baselineCats.append(baselineGalaxyCatalog(self.galaxyDBObj, obs_metadata=self.obs_metadata_compound,
                                                  column_outputs=['lsst_g', 'lsst_i', 'sigma_lsst_g', 'sigma_lsst_i']))

        testName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                'testSetup_testSetupCatUncertainty.txt')

        baseName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'testSetup_baseSetupCatUncertainty.txt')

        if os.path.exists(testName):
            os.unlink(testName)

        if os.path.exists(baseName):
            os.unlink(baseName)

        basedtype = np.dtype([('raObserved', np.float), ('decObserved', np.float),
                                 ('lsst_g', np.float), ('lsst_i', np.float),
                                 ('sigma_lsst_g',np.float), ('sigma_lsst_i', np.float)])

        for (testCatClass, dbo, baselineCat, msgr) in zip(testCatClasses, testCatDBs, baselineCats, msgroot):

            testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata,
                                             dbConnection=dbo,
                                             catalogClass=testCatClass,
                                             uncertainty=True)

            testdtype = np.dtype([('raObserved', np.float), ('decObserved', np.float),
                                     ('lsst_g', np.float), ('sigma_lsst_g', np.float)])

            testCat.write_catalog(testName)
            baselineCat.write_catalog(baseName)

            testData = np.genfromtxt(testName, dtype=testdtype, delimiter=',')
            baseData = np.genfromtxt(baseName, dtype=basedtype, delimiter=',')
            self.assertGreater(len(testData), 0)
            self.assertGreater(len(baseData), 0)

            ct = 0
            for b, t in zip(baseData, testData):
                self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12,
                                       msg = '%s single column; %.12e != %.12e ' % (msgr, b['lsst_g'], t['lsst_g']))
                self.assertAlmostEqual(b['sigma_lsst_g'], t['sigma_lsst_g'], 12,
                                       msg = '%s sigle column; %.12e != %.12e ' % (msgr, b['sigma_lsst_i'], t['sigma_lsst_g']))
                ct +=1

            self.assertGreater(ct, 0)

            testdtype = np.dtype([('raObserved', np.float), ('decObserved', np.float),
                                     ('lsst_g', np.float), ('sigma_lsst_g', np.float),
                                     ('lsst_i', np.float), ('sigma_lsst_i', np.float)])

            testCat = setupPhotometryCatalog(obs_metadata=self.obs_metadata_compound,
                                             dbConnection=dbo,
                                             catalogClass=testCatClass,
                                             uncertainty=True)
            testCat.write_catalog(testName)
            testData = np.genfromtxt(testName, dtype=testdtype, delimiter=',')
            self.assertGreater(len(testData), 0)
            ct = 0
            for b, t in zip(baseData, testData):
                self.assertAlmostEqual(b['lsst_g'], t['lsst_g'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['lsst_g'], t['lsst_g']))
                self.assertAlmostEqual(b['lsst_i'], t['lsst_i'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['lsst_i'], t['lsst_i']))
                self.assertAlmostEqual(b['sigma_lsst_g'], t['sigma_lsst_g'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['sigma_lsst_g'], t['lsst_g']))
                self.assertAlmostEqual(b['sigma_lsst_i'], t['sigma_lsst_i'], 12,
                                       msg = '%s double column; %.12e != %.12e ' % (msgr, b['sigma_lsst_i'], t['sigma_lsst_i']))
                ct +=1

            self.assertGreater(ct, 0)

            if os.path.exists(testName):
                os.unlink(testName)
            if os.path.exists(baseName):
                os.unlink(baseName)

class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
