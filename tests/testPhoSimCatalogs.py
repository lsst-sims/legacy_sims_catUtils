from __future__ import with_statement
import os
import numpy as np
import unittest
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.utils import defaultSpecMap, altAzPaFromRaDec, ObservationMetaData
from lsst.sims.catalogs.definitions import CompoundInstanceCatalog
from lsst.sims.catalogs.db import fileDBObject
from lsst.sims.catUtils.baseCatalogModels import SNDBObj
from lsst.sims.catUtils.utils import (testStarsDBObj, testGalaxyDiskDBObj,
                                      testGalaxyBulgeDBObj, testGalaxyAgnDBObj)
from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogSersic2D, PhoSimCatalogPoint,
                                                          PhoSimCatalogZPoint, PhoSimCatalogSSM,
                                                          PhoSimCatalogSN)
from lsst.sims.catUtils.utils import makePhoSimTestDB


def createTestSNDB():
    """
    Create a CatalogDBObject-like database that contains everything needed to create
    a supernova catalog.  Return the CatalogDBObject and an ObservationMetaData pointing
    to its center.

    Note: the OpsimMetaData for the returned ObservationMetaData will be inconsistent.
    It is just there so that PhoSim InstanceCatalog classes can write out a header.
    """

    raCenter = 23.0
    decCenter = -19.0
    radius = 0.1
    obs = ObservationMetaData(pointingRA=raCenter, pointingDec=decCenter, boundType='circle',
                              boundLength=radius, rotSkyPos=33.0, mjd=59580.0,
                              bandpassName='r')

    # these will be totally inconsistent; just need something to put in the header
    obs.OpsimMetaData = {'dist2Moon': 0.1, 'moonalt': -0.2,
                         'moonra': 1.1, 'moondec': 0.5,
                         'rottelpos': 0.4, 'sunalt': -1.1}

    rng = np.random.RandomState(88)
    n_obj = 10
    rr = rng.random_sample(n_obj)*radius
    theta = rng.random_sample(n_obj)*2.0*np.pi
    ra_list = raCenter + rr*np.cos(theta)
    dec_list = decCenter + rr*np.sin(theta)
    t0_list = 0.5*rng.random_sample(n_obj)
    x0_list = rng.random_sample(n_obj)*1.0e-4
    x1_list = rng.random_sample(n_obj)
    c_list = rng.random_sample(n_obj)*2.0
    z_list = rng.random_sample(n_obj)*1.1 + 0.1

    txt_file_name = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                 'test_phosim_sn_source.txt')

    if os.path.exists(txt_file_name):
        os.unlink(txt_file_name)

    dtype = np.dtype([('id', int), ('snra', float), ('sndec', float),
                      ('t0', float), ('x0', float), ('x1', float),
                      ('c', float), ('redshift', float), ('galtileid', int)])

    with open(txt_file_name, 'w') as output_file:
        for ix, (ra, dec, t0, x0, x1, c, z) in enumerate(zip(ra_list, dec_list, t0_list, x0_list,
                                                             x1_list, c_list, z_list)):
            output_file.write('%d %e %e %e %e %e %e %e %d\n' %
                              (ix, ra, dec, t0, x0, x1, c, z, ix))

    class testSNDBObj(SNDBObj, fileDBObject):
        tableid = 'test'

        def query_columns(self, *args, **kwargs):
            return fileDBObject.query_columns(self, *args, **kwargs)

    dbobj = testSNDBObj(txt_file_name, runtable='test', dtype=dtype)

    if os.path.exists(txt_file_name):
        os.unlink(txt_file_name)

    return dbobj, obs


test_header_map = {'dist2moon': ('dist2moon', np.degrees),
                   'moonalt': ('moonalt', np.degrees),
                   'moondec': ('moondec', np.degrees),
                   'moonra': ('moonra', np.degrees),
                   'rottelpos': ('rottelpos', np.degrees),
                   'sunalt': ('sunalt', np.degrees)}


def setup_module(module):
    lsst.utils.tests.init()


class PhoSimCatalogTest(unittest.TestCase):

    longMessage = True

    def setUp(self):
        self.obs_metadata = makePhoSimTestDB(size=10)
        self.bulgeDB = testGalaxyBulgeDBObj(driver='sqlite', database='PhoSimTestDatabase.db')
        self.diskDB = testGalaxyDiskDBObj(driver='sqlite', database='PhoSimTestDatabase.db')
        self.agnDB = testGalaxyAgnDBObj(driver='sqlite', database='PhoSimTestDatabase.db')
        self.starDB = testStarsDBObj(driver='sqlite', database='PhoSimTestDatabase.db')
        filter_translation = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4, 'y': 5}
        alt, az, pa = altAzPaFromRaDec(self.obs_metadata.pointingRA,
                                       self.obs_metadata.pointingDec,
                                       self.obs_metadata, includeRefraction=False)
        self.control_header = ['moondec %.7f\n' % np.degrees(self.obs_metadata.OpsimMetaData['moondec']),
                               'rottelpos %.7f\n' % np.degrees(self.obs_metadata.OpsimMetaData['rottelpos']),
                               'declination %.7f\n' % self.obs_metadata.pointingDec,
                               'moonalt %.7f\n' % np.degrees(self.obs_metadata.OpsimMetaData['moonalt']),
                               'rotskypos %.7f\n' % self.obs_metadata.rotSkyPos,
                               'moonra %.7f\n' % np.degrees(self.obs_metadata.OpsimMetaData['moonra']),
                               'sunalt %.7f\n' % np.degrees(self.obs_metadata.OpsimMetaData['sunalt']),
                               'mjd %.7f\n' % (self.obs_metadata.mjd.TAI+16.5/86400.0),
                               'azimuth %.7f\n' % az,
                               'rightascension %.7f\n' % self.obs_metadata.pointingRA,
                               'dist2moon %.7f\n' % np.degrees(self.obs_metadata.OpsimMetaData['dist2moon']),
                               'filter %d\n' % filter_translation[self.obs_metadata.bandpass],
                               'altitude %.7f\n' % alt]

    def tearDown(self):
        del self.starDB
        del self.bulgeDB
        del self.diskDB
        del self.agnDB
        if os.path.exists('PhoSimTestDatabase.db'):
            os.unlink('PhoSimTestDatabase.db')

    def verify_catalog(self, file_name):
        """
        Verify that the catalog specified by file_name has the expected header and is not empty
        """
        with open(file_name, 'r') as input_file:
            lines = input_file.readlines()
            for control_line in self.control_header:
                self.assertIn(control_line, lines)

            self.assertGreater(len(lines), len(self.control_header)+3)

    def testSpecFileMap(self):
        """
        Test that the PhoSim InstanceCatalog SpecFileMaps map MLT dwarf spectra
        to the starSED/phoSimMLT/ directory (that is where the MLT spectra which
        have been clipped to honor PhoSim's 'no more than 24,999 lines per SED
        file' requirement reside)
        """

        cat = PhoSimCatalogPoint(self.starDB, obs_metadata=self.obs_metadata)
        self.assertEqual(cat.specFileMap['lte_123.txt'], 'starSED/phoSimMLT/lte_123.txt.gz')
        self.assertEqual(cat.specFileMap['lte_123.txt.gz'], 'starSED/phoSimMLT/lte_123.txt.gz')
        self.assertNotEqual(cat.specFileMap['lte_123.txt'], defaultSpecMap['lte_123.txt'])

        # verify that the usual stellar mappings still apply
        for test_name in ('kp_123.txt', 'km_123.txt', 'Const_123.txt', 'Exp_123.txt', 'Burst_123.txt',
                          'bergeron_123.txt', 'burrows_123.txt', 'm1_123.txt', 'L1_123.txt', 'l1_123.txt',
                          'Inst_123.txt'):

            self.assertEqual(cat.specFileMap[test_name], defaultSpecMap[test_name])

    def test_incomplete_obs(self):
        """
        Test that an exception gets raised if you try to make a PhoSim InstanceCatalog
        with an ObservationMetaData that lacks RA, Dec, mjd, bandpass, or rotSkyPos
        """
        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'bad_obs_test_phosim_cat.txt')
        obs = ObservationMetaData(pointingDec=19.0, mjd=43000.0, rotSkyPos=19.0, bandpassName='u')
        with self.assertRaises(TypeError):
            cat = PhoSimCatalogPoint(self.starDB, obs_metadata=obs)
            cat.phoSimHeaderMap = {}
            cat.write_catalog(catName)

        obs = ObservationMetaData(pointingRA=19.0, mjd=43000.0, rotSkyPos=19.0, bandpassName='u')
        with self.assertRaises(TypeError):
            cat = PhoSimCatalogPoint(self.starDB, obs_metadata=obs)
            cat.phoSimHeaderMap = {}
            cat.write_catalog(catName)

        obs = ObservationMetaData(pointingRA=88.0, pointingDec=19.0, rotSkyPos=19.0, bandpassName='u')
        with self.assertRaises(RuntimeError):
            cat = PhoSimCatalogPoint(self.starDB, obs_metadata=obs)
            cat.phoSimHeaderMap = {}
            cat.write_catalog(catName)

        obs = ObservationMetaData(pointingRA=88.0, pointingDec=19.0, mjd=43000.0, bandpassName='u')
        with self.assertRaises(TypeError):
            cat = PhoSimCatalogPoint(self.starDB, obs_metadata=obs)
            cat.phoSimHeaderMap = {}
            cat.write_catalog(catName)

        obs = ObservationMetaData(pointingRA=88.0, pointingDec=19.0, mjd=43000.0, rotSkyPos=19.0)
        with self.assertRaises(KeyError):
            cat = PhoSimCatalogPoint(self.starDB, obs_metadata=obs)
            cat.phoSimHeaderMap = {}
            cat.write_catalog(catName)

        obs = ObservationMetaData(pointingRA=88.0, pointingDec=19.0, mjd=43000.0, rotSkyPos=19.0,
                                  bandpassName='u')

        cat = PhoSimCatalogPoint(self.starDB, obs_metadata=obs)
        cat.phoSimHeaderMap = {}
        cat.write_catalog(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testCatalog(self):
        """
        This test writes a PhoSim input catalog and compares it, one line at a time
        to a previously written catalog that should be identical.
        """
        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata = self.obs_metadata)
        testDisk = PhoSimCatalogSersic2D(self.diskDB, obs_metadata = self.obs_metadata)
        testAgn = PhoSimCatalogZPoint(self.agnDB, obs_metadata = self.obs_metadata)
        testStar = PhoSimCatalogPoint(self.starDB, obs_metadata = self.obs_metadata)

        catName = os.path.join(getPackageDir('sims_catUtils'),
                               'tests', 'scratchSpace', 'phoSimTestCatalog.txt')

        testBulge.phoSimHeaderMap = test_header_map
        testBulge.write_catalog(catName)
        testDisk.write_catalog(catName, write_header=False, write_mode='a')
        testAgn.write_catalog(catName, write_header=False, write_mode='a')
        testStar.write_catalog(catName, write_header=False, write_mode='a')

        self.verify_catalog(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testHeaderMap(self):
        """
        Test the behavior of the phoSimHeaderMap
        """
        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=self.obs_metadata)
        testBulge.phoSimHeaderMap = {'lunar_distance': ('dist2moon', None),
                                     'rotation_of_the_telescope': ('rottelpos', np.degrees),
                                     'other_rotation': ('rottelpos', lambda x: x*x)}

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'header_map_phosim_catalog.txt')
        testBulge.write_catalog(catName)

        with open(catName, 'r') as input_file:
            input_header = {}
            for line in input_file:
                vv = line.split()
                if vv[0] != 'object':
                    input_header[vv[0]] = vv[1]
                else:
                    break

        self.assertIn('rightascension', input_header)
        self.assertIn('declination', input_header)
        self.assertIn('altitude', input_header)
        self.assertIn('azimuth', input_header)
        self.assertIn('filter', input_header)
        self.assertIn('rotskypos', input_header)
        self.assertIn('mjd', input_header)
        self.assertIn('lunar_distance', input_header)
        self.assertAlmostEqual(float(input_header['lunar_distance']),
                               self.obs_metadata.OpsimMetaData['dist2moon'], 6)
        self.assertIn('rotation_of_the_telescope', input_header)
        self.assertAlmostEqual(float(input_header['rotation_of_the_telescope']),
                               np.degrees(self.obs_metadata.OpsimMetaData['rottelpos']),
                               delta=1.0e-6*np.degrees(self.obs_metadata.OpsimMetaData['rottelpos']))
        self.assertIn('other_rotation', input_header)
        self.assertAlmostEqual(float(input_header['other_rotation']),
                               self.obs_metadata.OpsimMetaData['rottelpos']**2,
                               delta=1.0e-6*self.obs_metadata.OpsimMetaData['rottelpos']**2)
        self.assertEqual(len(input_header), 10)

        if os.path.exists(catName):
            os.unlink(catName)

    def testBlankHeaderMap(self):
        """
        Test behavior of a blank header map
        """
        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=self.obs_metadata)
        testBulge.phoSimHeaderMap = {}

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'blank_header_map_phosim_catalog.txt')
        testBulge.write_catalog(catName)

        with open(catName, 'r') as input_file:
            input_header = {}
            for line in input_file:
                vv = line.split()
                if vv[0] != 'object':
                    input_header[vv[0]] = vv[1]
                else:
                    break

        # verify that only the default header parameters are included in the
        # PhoSimInstanceCatalog, even though obs_metadata has a non-None
        # OpsimMetaData
        self.assertIn('rightascension', input_header)
        self.assertIn('declination', input_header)
        self.assertIn('altitude', input_header)
        self.assertIn('azimuth', input_header)
        self.assertIn('filter', input_header)
        self.assertIn('rotskypos', input_header)
        self.assertIn('mjd', input_header)
        self.assertEqual(len(input_header), 7)
        self.assertGreater(len(self.obs_metadata.OpsimMetaData), 0)

        if os.path.exists(catName):
            os.unlink(catName)

    def testNoHeaderMap(self):
        """
        Test that the correct error is raised if no header map is specified
        """
        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=self.obs_metadata)

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'no_header_map_phosim_catalog.txt')

        with self.assertRaises(RuntimeError) as context:
            testBulge.write_catalog(catName)

        self.assertIn("without specifying a phoSimHeaderMap",
                      context.exception.args[0])

        # now make sure that the exception is raised, even if ObservationMetaData
        # does not have an OpsimMetaData
        obs = ObservationMetaData(pointingRA=35.0, pointingDec=-23.0,
                                  mjd=43900.0, rotSkyPos=22.0,
                                  boundType='circle', boundLength=1.75)

        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=obs)
        with self.assertRaises(RuntimeError) as context:
            testBulge.write_catalog(catName)

        self.assertIn("without specifying a phoSimHeaderMap",
                      context.exception.args[0])
        self.assertIn("you may wish to consider adding default PhoSim parameters",
                      context.exception.args[0])

        if os.path.exists(catName):
            os.unlink(catName)

    def test_default_values_in_header_map(self):
        """
        Test that default PhoSim header values in the header map get appropriately applied
        """
        test_header_map = {'lunar_distance': ('dist2moon', None),
                           'nsnap': 3}

        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=self.obs_metadata)
        testBulge.phoSimHeaderMap = test_header_map

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'default_value_header_map_phosim_catalog.txt')
        testBulge.write_catalog(catName)

        with open(catName, 'r') as input_file:
            input_header = {}
            for line in input_file:
                vv = line.split()
                if vv[0] != 'object':
                    input_header[vv[0]] = vv[1]
                else:
                    break

        self.assertIn('rightascension', input_header)
        self.assertIn('declination', input_header)
        self.assertIn('altitude', input_header)
        self.assertIn('azimuth', input_header)
        self.assertIn('filter', input_header)
        self.assertIn('rotskypos', input_header)
        self.assertIn('mjd', input_header)
        self.assertIn('lunar_distance', input_header)
        self.assertAlmostEqual(float(input_header['lunar_distance']),
                               self.obs_metadata.OpsimMetaData['dist2moon'], 6)
        self.assertIn('nsnap', input_header)
        self.assertEqual(int(input_header['nsnap']), 3)
        self.assertEqual(len(input_header), 9)

        if os.path.exists(catName):
            os.unlink(catName)

    def test_non_existent_values_in_header_map(self):
        """
        Test that header params that are defined in the header map but not
        in OpsimMetaData are ommitted from the header
        """
        test_header_map = {'lunar_distance': ('dist2moon', None),
                           'nsnap': 3,
                           'nonesense': ('gobbledygook', lambda x: 2.0*x)}

        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=self.obs_metadata)
        testBulge.phoSimHeaderMap = test_header_map

        catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'nonexistent_value_header_map_phosim_catalog.txt')
        testBulge.write_catalog(catName)

        with open(catName, 'r') as input_file:
            input_header = {}
            for line in input_file:
                vv = line.split()
                if vv[0] != 'object':
                    input_header[vv[0]] = vv[1]
                else:
                    break

        self.assertIn('rightascension', input_header)
        self.assertIn('declination', input_header)
        self.assertIn('altitude', input_header)
        self.assertIn('azimuth', input_header)
        self.assertIn('filter', input_header)
        self.assertIn('rotskypos', input_header)
        self.assertIn('mjd', input_header)
        self.assertIn('lunar_distance', input_header)
        self.assertAlmostEqual(float(input_header['lunar_distance']),
                               self.obs_metadata.OpsimMetaData['dist2moon'], 6)
        self.assertIn('nsnap', input_header)
        self.assertEqual(int(input_header['nsnap']), 3)
        self.assertEqual(len(input_header), 9)

        # now try it with no OpsimMetaData at all
        obs = ObservationMetaData(pointingRA=23.0, pointingDec=-11.0,
                                  mjd=43000.0, rotSkyPos=44.0,
                                  bandpassName='g',
                                  boundType='circle', boundLength=1.0)
        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=obs)
        testBulge.phoSimHeaderMap = test_header_map
        testBulge.write_catalog(catName)

        with open(catName, 'r') as input_file:
            input_header = {}
            for line in input_file:
                vv = line.split()
                if vv[0] != 'object':
                    input_header[vv[0]] = vv[1]
                else:
                    break

        self.assertIn('rightascension', input_header)
        self.assertIn('declination', input_header)
        self.assertIn('altitude', input_header)
        self.assertIn('azimuth', input_header)
        self.assertIn('filter', input_header)
        self.assertIn('rotskypos', input_header)
        self.assertIn('mjd', input_header)
        self.assertIn('nsnap', input_header)
        self.assertEqual(int(input_header['nsnap']), 3)
        self.assertEqual(len(input_header), 8)

        if os.path.exists(catName):
            os.unlink(catName)

    def testCompoundCatalog(self):
        """
        This test writes a PhoSim input catalog and compares it, one line at a time
        to a previously written catalog that should be identical.

        This test uses CompoundInstanceCatalog
        """

        # first, generate the catalog without a CompoundInstanceCatalog
        single_catName = os.path.join(getPackageDir('sims_catUtils'),
                                      'tests', 'scratchSpace', 'phoSimTestCatalog_single.txt')

        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata = self.obs_metadata)
        testDisk = PhoSimCatalogSersic2D(self.diskDB, obs_metadata = self.obs_metadata)
        testAgn = PhoSimCatalogZPoint(self.agnDB, obs_metadata = self.obs_metadata)
        testStar = PhoSimCatalogPoint(self.starDB, obs_metadata = self.obs_metadata)

        testBulge.phoSimHeaderMap = test_header_map
        testBulge.write_catalog(single_catName)
        testDisk.write_catalog(single_catName, write_header=False, write_mode='a')
        testAgn.write_catalog(single_catName, write_header=False, write_mode='a')
        testStar.write_catalog(single_catName, write_header=False, write_mode='a')

        # now, generate the catalog using CompoundInstanceCatalog
        #
        # because the CompoundCatalogDBObject requires that database
        # connection parameters be set in the input CatalogDBObject
        # daughter class definitions, we have to declare dummy
        # CatalogDBObject daughter classes below

        class dummyDBbase(object):
            driver = 'sqlite'
            database = 'PhoSimTestDatabase.db'

        class dummyBulgeDB(dummyDBbase, testGalaxyBulgeDBObj):
            objid = 'dummy_bulge'

        class dummyDiskDB(dummyDBbase, testGalaxyDiskDBObj):
            objid = 'dummy_disk'

        class dummyAgnDB(dummyDBbase, testGalaxyAgnDBObj):
            objid = 'dummy_agn'

        class dummyStarDB(dummyDBbase, testStarsDBObj):
            objid = 'dummy_stars'

        compoundCatalog = CompoundInstanceCatalog([PhoSimCatalogSersic2D, PhoSimCatalogSersic2D,
                                                   PhoSimCatalogZPoint, PhoSimCatalogPoint],
                                                  [dummyBulgeDB, dummyDiskDB, dummyAgnDB, dummyStarDB],
                                                  obs_metadata=self.obs_metadata)

        self.assertEqual(len(compoundCatalog._dbObjectGroupList[0]), 3)

        compound_catName = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                        'phoSimTestCatalog_compound.txt')

        compoundCatalog.phoSimHeaderMap = test_header_map
        compoundCatalog.write_catalog(compound_catName)

        # verify that the compound catalog is what we expect
        self.verify_catalog(compound_catName)

        # verify that the two catalogs are equivalent
        with open(single_catName, 'r') as single_file:
            with open(compound_catName, 'r') as compound_file:
                single_lines = single_file.readlines()
                compound_lines = compound_file.readlines()

                for line in single_lines:
                    self.assertIn(line, compound_lines)

                for line in compound_lines:
                    self.assertIn(line, single_lines)

        if os.path.exists(compound_catName):
            os.unlink(compound_catName)

        if os.path.exists(single_catName):
            os.unlink(single_catName)

    def testPointSourceSchema(self):
        """
        Create a PhoSim InstanceCatalog of point sources (stars).  Verify
        that the schema of the actual objects conforms to what PhoSim expects,
        as defined here

        https://bitbucket.org/phosim/phosim_release/wiki/Instance%20Catalog
        """
        cat = PhoSimCatalogPoint(self.starDB, obs_metadata=self.obs_metadata)
        cat.phoSimHeaderMap = test_header_map
        cat_name = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                'phosim_point_source_schema_cat.txt')
        if os.path.exists(cat_name):
            os.unlink(cat_name)

        cat.write_catalog(cat_name)

        with open(cat_name, 'r') as input_file:
            cat_lines = input_file.readlines()

        n_obj = 0
        for line in cat_lines:
            params = line.split()
            if len(params) > 2:
                n_obj += 1
                self.assertEqual(len(params), 17)
                self.assertEqual(params[0], 'object')
                self.assertEqual(round(float(params[1])), float(params[1]), 10)  # id
                float(params[2])  # ra
                float(params[3])  # dec
                float(params[4])  # mag norm
                self.assertIn('starSED', params[5])  # sed name
                self.assertAlmostEqual(float(params[6]), 0.0, 10)  # redshift
                self.assertAlmostEqual(float(params[7]), 0.0, 10)  # gamma1
                self.assertAlmostEqual(float(params[8]), 0.0, 10)  # gamma2
                self.assertAlmostEqual(float(params[9]), 0.0, 10)  # kappa
                self.assertAlmostEqual(float(params[10]), 0.0, 10)  # delta_ra
                self.assertAlmostEqual(float(params[11]), 0.0, 10)  # delta_dec
                self.assertEqual(params[12], 'point')  # source type
                dust_msg = ('It is possible you are outputting Milky Way dust parameters before '
                            'internal dust parameters; internal dust should come first')
                self.assertEqual(params[13], 'none', msg=dust_msg)  # internal dust
                self.assertEqual(params[14], 'CCM', msg=dust_msg)  # Milky Way dust
                self.assertGreater(float(params[15]), 0.0, msg=dust_msg)  # Av
                self.assertAlmostEqual(float(params[16]), 3.1, msg=dust_msg)  # Rv

        self.assertGreater(n_obj, 0)

        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testSersicSchema(self):
        """
        Create a PhoSim InstanceCatalog of Sersic profiles (galaxy bulges).  Verify
        that the schema of the actual objects conforms to what PhoSim expects,
        as defined here

        https://bitbucket.org/phosim/phosim_release/wiki/Instance%20Catalog
        """
        cat = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=self.obs_metadata)
        cat.phoSimHeaderMap = test_header_map
        cat_name = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                'phosim_sersic_schema_cat.txt')
        if os.path.exists(cat_name):
            os.unlink(cat_name)

        cat.write_catalog(cat_name)

        with open(cat_name, 'r') as input_file:
            cat_lines = input_file.readlines()

        n_obj = 0
        for line in cat_lines:
            params = line.split()
            if len(params) > 2:
                n_obj += 1
                self.assertEqual(len(params), 23)
                self.assertEqual(params[0], 'object')
                self.assertEqual(round(float(params[1])), float(params[1]), 10)  # id
                float(params[2])  # ra
                float(params[3])  # dec
                float(params[4])  # mag norm
                self.assertIn('galaxySED', params[5])  # sed name
                self.assertGreater(float(params[6]), 0.0, 10)  # redshift
                self.assertAlmostEqual(float(params[7]), 0.0, 10)  # gamma1
                self.assertAlmostEqual(float(params[8]), 0.0, 10)  # gamma2
                self.assertAlmostEqual(float(params[9]), 0.0, 10)  # kappa
                self.assertAlmostEqual(float(params[10]), 0.0, 10)  # delta_ra
                self.assertAlmostEqual(float(params[11]), 0.0, 10)  # delta_dec
                self.assertEqual(params[12], 'sersic2d')  # source type
                self.assertGreater(float(params[13]), 0.0)  # major axis
                self.assertGreater(float(params[14]), 0.0)  # minor axis
                self.assertGreater(float(params[15]), 0.0)  # position angle
                self.assertAlmostEqual(float(params[16]), 4.0, 13)  # n_s (bulges have sersic index=4)
                self.assertEqual(params[17], 'CCM')  # internal dust
                dust_msg = ('It is possible you are outputting Milky Way dust parameters before '
                            'internal dust parameters; internal dust should come first')
                self.assertLess(float(params[18]), 0.31, msg=dust_msg)  # Av
                self.assertLess(float(params[19]), 2.11, msg=dust_msg)  # Rv
                self.assertEqual(params[20], 'CCM')  # Milky Way dust
                self.assertGreater(float(params[21]), 0.0, msg=dust_msg)  # Av
                self.assertAlmostEqual(float(params[22]), 3.1, msg=dust_msg)  # Rv

        self.assertGreater(n_obj, 0)

        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testZPointSourceSchema(self):
        """
        Create a PhoSim InstanceCatalog of extra-galactic point sources (agns).  Verify
        that the schema of the actual objects conforms to what PhoSim expects,
        as defined here

        https://bitbucket.org/phosim/phosim_release/wiki/Instance%20Catalog
        """
        cat = PhoSimCatalogZPoint(self.agnDB, obs_metadata=self.obs_metadata)
        cat.phoSimHeaderMap = test_header_map
        cat_name = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                'phosim_agn_schema_cat.txt')
        if os.path.exists(cat_name):
            os.unlink(cat_name)

        cat.write_catalog(cat_name)

        with open(cat_name, 'r') as input_file:
            cat_lines = input_file.readlines()

        n_obj = 0
        for line in cat_lines:
            params = line.split()
            if len(params) > 2:
                n_obj += 1
                self.assertEqual(len(params), 17)
                self.assertEqual(params[0], 'object')
                self.assertEqual(round(float(params[1])), float(params[1]), 10)  # id
                float(params[2])  # ra
                float(params[3])  # dec
                float(params[4])  # mag norm
                self.assertIn('agnSED', params[5])  # sed name
                self.assertGreater(float(params[6]), 0.0)  # redshift
                self.assertAlmostEqual(float(params[7]), 0.0, 10)  # gamma1
                self.assertAlmostEqual(float(params[8]), 0.0, 10)  # gamma2
                self.assertAlmostEqual(float(params[9]), 0.0, 10)  # kappa
                self.assertAlmostEqual(float(params[10]), 0.0, 10)  # delta_ra
                self.assertAlmostEqual(float(params[11]), 0.0, 10)  # delta_dec
                self.assertEqual(params[12], 'point')  # source type
                dust_msg = ('It is possible you are outputting Milky Way dust parameters before '
                            'internal dust parameters; internal dust should come first')
                self.assertEqual(params[13], 'none', msg=dust_msg)  # internal dust
                self.assertEqual(params[14], 'CCM', msg=dust_msg)  # Milky Way dust
                self.assertGreater(float(params[15]), 0.0, msg=dust_msg)  # Av
                self.assertAlmostEqual(float(params[16]), 3.1, msg=dust_msg)  # Rv

        self.assertGreater(n_obj, 0)

        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testSNSchema(self):
        """
        Create a PhoSim InstanceCatalog of supernovae.  Verify
        that the schema of the actual objects conforms to what PhoSim expects,
        as defined here

        https://bitbucket.org/phosim/phosim_release/wiki/Instance%20Catalog
        """
        db, obs = createTestSNDB()
        cat = PhoSimCatalogSN(db, obs_metadata=obs)
        cat.writeSEDFile = False
        cat.phoSimHeaderMap = test_header_map
        cat_name = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                'phosim_sne_schema_cat.txt')
        if os.path.exists(cat_name):
            os.unlink(cat_name)

        cat.write_catalog(cat_name)

        with open(cat_name, 'r') as input_file:
            cat_lines = input_file.readlines()

        n_obj = 0
        for line in cat_lines:
            params = line.split()
            if len(params) > 2:
                n_obj += 1
                self.assertEqual(len(params), 17)
                self.assertEqual(params[0], 'object')
                self.assertEqual(round(float(params[1])), float(params[1]), 10)  # id
                float(params[2])  # ra
                float(params[3])  # dec
                float(params[4])  # mag norm
                self.assertIn('%s.dat' % obs.bandpass, params[5])  # sed name
                self.assertGreater(float(params[6]), 0.0)  # redshift
                self.assertAlmostEqual(float(params[7]), 0.0, 10)  # gamma1
                self.assertAlmostEqual(float(params[8]), 0.0, 10)  # gamma2
                self.assertAlmostEqual(float(params[9]), 0.0, 10)  # kappa
                self.assertAlmostEqual(float(params[10]), 0.0, 10)  # delta_ra
                self.assertAlmostEqual(float(params[11]), 0.0, 10)  # delta_dec
                self.assertEqual(params[12], 'point')  # source type
                dust_msg = ('It is possible you are outputting Milky Way dust parameters before '
                            'internal dust parameters; internal dust should come first')
                self.assertEqual(params[13], 'none', msg=dust_msg)  # internal dust
                self.assertEqual(params[14], 'CCM', msg=dust_msg)  # Milky Way dust
                self.assertGreater(float(params[15]), 0.0, msg=dust_msg)  # Av
                self.assertAlmostEqual(float(params[16]), 3.1, msg=dust_msg)  # Rv

        self.assertGreater(n_obj, 0)

        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testSSMSchema(self):
        """
        Create a PhoSim InstanceCatalog of point sources (stars) formatted by the
        PhoSimCatalogSSM class.  Verify that the schema of the actual objects conforms
        to what PhoSim expects, as defined here

        https://bitbucket.org/phosim/phosim_release/wiki/Instance%20Catalog
        """
        cat = PhoSimCatalogSSM(self.starDB, obs_metadata=self.obs_metadata)
        cat.phoSimHeaderMap = test_header_map
        cat_name = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                'phosim_ssm_schema_cat.txt')
        if os.path.exists(cat_name):
            os.unlink(cat_name)

        cat.write_catalog(cat_name)

        with open(cat_name, 'r') as input_file:
            cat_lines = input_file.readlines()

        n_obj = 0
        for line in cat_lines:
            params = line.split()
            if len(params) > 2:
                n_obj += 1
                self.assertEqual(len(params), 15)
                self.assertEqual(params[0], 'object')
                self.assertEqual(round(float(params[1])), float(params[1]), 10)  # id
                float(params[2])  # ra
                float(params[3])  # dec
                float(params[4])  # mag norm
                self.assertIn('starSED', params[5])  # sed name
                self.assertAlmostEqual(float(params[6]), 0.0, 10)  # redshift
                self.assertAlmostEqual(float(params[7]), 0.0, 10)  # gamma1
                self.assertAlmostEqual(float(params[8]), 0.0, 10)  # gamma2
                self.assertAlmostEqual(float(params[9]), 0.0, 10)  # kappa
                self.assertAlmostEqual(float(params[10]), 0.0, 10)  # delta_ra
                self.assertAlmostEqual(float(params[11]), 0.0, 10)  # delta_dec
                self.assertEqual(params[12], 'point')  # source type
                self.assertEqual(params[13], 'none')  # internal dust
                self.assertEqual(params[14], 'none')  # Milky Way dust

        self.assertGreater(n_obj, 0)

        if os.path.exists(cat_name):
            os.unlink(cat_name)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
