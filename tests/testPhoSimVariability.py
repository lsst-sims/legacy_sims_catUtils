from builtins import zip
import os
import unittest
import lsst.utils.tests
import numpy as np
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.utils import (testStarsDBObj, testGalaxyDiskDBObj,
                                      testGalaxyBulgeDBObj, testGalaxyAgnDBObj)
from lsst.sims.catUtils.exampleCatalogDefinitions import (PhoSimCatalogSersic2D, PhoSimCatalogPoint,
                                                          PhoSimCatalogZPoint)
from lsst.sims.catUtils.utils import makePhoSimTestDB
from lsst.sims.catUtils.mixins import VariabilityStars, VariabilityGalaxies
from lsst.sims.catUtils.mixins import VariabilityAGN
from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels
from lsst.sims.catUtils.utils import TestVariabilityMixin
from lsst.sims.catUtils.mixins import AstrometryStars, AstrometryGalaxies


def setup_module(module):
    lsst.utils.tests.init()


class PhoSimPointVariable(PhoSimCatalogPoint, VariabilityStars, TestVariabilityMixin):
    catalog_type = __file__ + 'pho_sim_point_variable'
    pass


class PhoSimZPointVariable(PhoSimCatalogZPoint, VariabilityAGN, TestVariabilityMixin):
    catalog_type = __file__ + 'pho_sim_z_point_variable'
    pass


class AgnControlCatalog(InstanceCatalog, VariabilityGalaxies, TestVariabilityMixin, AstrometryGalaxies):
    catalog_type = __file__ + "agn_control_catalog"
    column_outputs = ['magNorm', 'delta_rAgn']


class BulgeControlCatalog(InstanceCatalog, AstrometryGalaxies):
    catalog_type = __file__ + "bulge_control_catalog"
    column_outputs = ['magNorm']


class DiskControlCatalog(InstanceCatalog, AstrometryGalaxies):
    catalog_type = __file__ + "disk_control_catalog"
    column_outputs = ['magNorm']


class StarControlCatalog(InstanceCatalog, AstrometryStars, VariabilityStars, TestVariabilityMixin):
    catalog_type = __file__ + "star_control_catalog"
    column_outputs = ['magNorm', 'delta_lsst_r']


class PhoSimVariabilityTest(unittest.TestCase):
    """
    This class will test that variability gets correctly propagated into
    PhoSim catalogs
    """

    longMessage = True

    @classmethod
    def setUpClass(cls):
        cls.dbName = 'PhoSimVariabilityDatabase.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        cls.obs_metadata = makePhoSimTestDB(size=10, filename=cls.dbName)

        cls.bulgeDB = testGalaxyBulgeDBObj(driver='sqlite', database=cls.dbName)
        cls.diskDB = testGalaxyDiskDBObj(driver='sqlite', database=cls.dbName)
        cls.agnDB = testGalaxyAgnDBObj(driver='sqlite', database=cls.dbName)
        cls.starDB = testStarsDBObj(driver='sqlite', database=cls.dbName)

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        del cls.bulgeDB
        del cls.diskDB
        del cls.agnDB
        del cls.starDB
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

    def testAgn(self):
        """
        Test that variability is correctly added to PhoSim Agn catalogs
        by outputting both a variable PhoSim catalog and a control catalog
        and making sure that the magNorm column in the PhoSim catalog
        is equal to the sum of the magNorm column in the control plus
        the detla_mag column from Variability.
        """
        baseline = AgnControlCatalog(self.agnDB, obs_metadata=self.obs_metadata)
        test = PhoSimZPointVariable(self.agnDB, obs_metadata=self.obs_metadata)

        for bb, tt in zip(baseline.iter_catalog(), test.iter_catalog()):
            msg = 'baseline mag %.6e; delta %.6e' % (bb[0], bb[1])
            self.assertAlmostEqual(bb[0] + bb[1], tt[4], 10, msg=msg)
            self.assertGreater(np.abs(bb[1]), 0.0)

    def testStars(self):
        """
        Test that variability is correctly added to PhoSim star catalogs
        by outputting both a variable PhoSim catalog and a control catalog
        and making sure that the magNorm column in the PhoSim catalog
        is equal to the sum of the magNorm column in the control plus
        the detla_mag column from Variability.
        """
        baseline = StarControlCatalog(self.starDB, obs_metadata=self.obs_metadata)
        test = PhoSimPointVariable(self.starDB, obs_metadata=self.obs_metadata)

        for bb, tt in zip(baseline.iter_catalog(), test.iter_catalog()):
            msg = 'baseline mag %.6e; delta %.6e' % (bb[0], bb[1])
            self.assertAlmostEqual(bb[0] + bb[1], tt[4], 10, msg=msg)
            self.assertGreater(np.abs(bb[1]), 0.0)

    def testBulges(self):
        """
        Make sure that the magNorm output to PhoSim catalogs that lack
        variability is the same as the column 'magNorm' taken from the database
        """
        baseline = BulgeControlCatalog(self.bulgeDB, obs_metadata=self.obs_metadata)
        test = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=self.obs_metadata)

        for bb, tt in zip(baseline.iter_catalog(), test.iter_catalog()):
            self.assertAlmostEqual(bb[0], tt[4], 10)

    def testDisks(self):
        baseline = DiskControlCatalog(self.diskDB, obs_metadata=self.obs_metadata)
        test = PhoSimCatalogSersic2D(self.diskDB, obs_metadata=self.obs_metadata)

        for bb, tt in zip(baseline.iter_catalog(), test.iter_catalog()):
            self.assertAlmostEqual(bb[0], tt[4], 10)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
