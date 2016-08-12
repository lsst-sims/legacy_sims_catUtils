from __future__ import with_statement
import os
import numpy
import unittest
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.utils import makeStarTestDB, myTestStars
from lsst.sims.catalogs.utils import makeGalTestDB, myTestGals
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import compound
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies
from lsst.sims.photUtils import BandpassDict

class FakeStellarVariabilityMixin(object):

    @compound('delta_test_u', 'delta_test_g', 'delta_test_r')
    def get_variability(self):
        ra = self.column_by_name('raJ2000')
        return numpy.array([4.0*numpy.ones(len(ra)),
                            10.0*numpy.ones(len(ra)),
                            20.0*numpy.ones(len(ra))])


class StellarBaselineCatalogClass(InstanceCatalog, PhotometryStars):

    default_columns = [('galacticAv', 0.1, float)]

    def get_sedFilename(self):

        star_seds = ['km20_5750.fits_g40_5790',
                     'kp10_9250.fits_g40_9250',
                     'bergeron_6500_85.dat_6700']

        ra = self.column_by_name('raJ2000')
        return numpy.array([star_seds[i%3] for i in range(len(ra))])

    @compound('test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y')
    def get_test_mags(self):
        if not hasattr(self, 'variabilitybandpassDict'):
            self.variabilityBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        self._loadSedList(self.variabilityBandpassDict.wavelenMatch)
        if not hasattr(self, '_sedList'):
            return numpy.ones((6,0))

        return self._magnitudeGetter(self.variabilityBandpassDict, self.get_test_mags._colnames)


class StellarVariabilityCatalogClass(StellarBaselineCatalogClass, FakeStellarVariabilityMixin):
    pass


class FakeGalaxyVariabilityMixin(object):

    @compound('delta_test_agn_u', 'delta_test_agn_g', 'delta_test_agn_r',
              'delta_test_bulge_u', 'delta_test_disk_i')
    def get_variability(self):
        ra = self.column_by_name('raJ2000')
        return numpy.array([1.0*numpy.ones(len(ra)),
                            2.0*numpy.ones(len(ra)),
                            3.0*numpy.ones(len(ra)),
                            4.0*numpy.ones(len(ra)),
                            5.0*numpy.ones(len(ra))])


class GalaxyBaselineCatalogClass(InstanceCatalog, PhotometryGalaxies):

    @compound('internalAvBulge', 'internalAvDisk')
    def get_internalAv(self):
        ra = self.column_by_name('raJ2000')
        return numpy.array([2.5*numpy.ones(len(ra)), 2.5*numpy.ones(len(ra))])

    @compound('sedFilenameBulge', 'sedFilenameDisk', 'sedFilenameAgn')
    def get_filenames(self):
        ra = self.column_by_name('raJ2000')

        galaxy_seds = ['Const.80E07.02Z.spec','Inst.80E07.002Z.spec','Burst.19E07.0005Z.spec']
        agn_sed = 'agn.spec'

        agnSeds = []
        for ii in range(len(ra)):
            agnSeds.append(agn_sed)

        bulgeSeds = [galaxy_seds[(ii+1)%3] for ii in range(len(ra))]
        diskSeds = [galaxy_seds[ii%3] for ii in range(len(ra))]

        return numpy.array([bulgeSeds, diskSeds, agnSeds])

    @compound('test_bulge_u', 'test_bulge_g', 'test_bulge_r',
              'test_bulge_i' ,'test_bulge_z', 'test_bulge_y')
    def get_test_bulge_mags(self):

        if not hasattr(self, 'testBandpassDict'):
            self.testBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        return self._magnitudeGetter('bulge', self.testBandpassDict,
                                     self.get_test_bulge_mags._colnames)


    @compound('test_disk_u', 'test_disk_g', 'test_disk_r',
              'test_disk_i', 'test_disk_z', 'test_disk_y')
    def get_test_disk_mags(self):

        if not hasattr(self, 'testBandpassDict'):
            self.testBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        return self._magnitudeGetter('disk', self.testBandpassDict,
                                     self.get_test_disk_mags._colnames)


    @compound('test_agn_u', 'test_agn_g', 'test_agn_r',
              'test_agn_i', 'test_agn_z', 'test_agn_y')
    def get_test_agn_mags(self):

        if not hasattr(self, 'testBandpassDict'):
            self.testBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        return self._magnitudeGetter('agn', self.testBandpassDict,
                                     self.get_test_agn_mags._colnames)


    @compound('test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y')
    def get_test_total_mags(self):
        idList = self.column_by_name('uniqueId')
        numObj = len(idList)
        output = []
        for columnName in self.get_test_total_mags._colnames:
            if columnName not in self._actually_calculated_columns:
                sub_list = [numpy.NaN]*numObj
            else:
                bandpass = columnName[-1]
                bulge = self.column_by_name('test_bulge_%s' % bandpass)
                disk = self.column_by_name('test_disk_%s' % bandpass)
                agn = self.column_by_name('test_agn_%s' % bandpass)
                sub_list = self.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

            output.append(sub_list)
        return numpy.array(output)


class GalaxyVariabilityCatalogClass(GalaxyBaselineCatalogClass, FakeGalaxyVariabilityMixin):
    pass

class VariabilityDesignTest(unittest.TestCase):
    """
    This unit test case will test that the general
    variability design was correclty implemented.
    It will not test an particular model of variability.
    It will merely test that, given a mean and delta magnitude,
    the InstanceCatalog class can correctly calculate the varying
    magnitude of an object.
    """

    @classmethod
    def setUpClass(cls):
        cls.starDbName = 'VariabilityInfrastructureTestStarDB.db'
        if os.path.exists(cls.starDbName):
            os.unlink(cls.starDbName)

        makeStarTestDB(filename=cls.starDbName)

        cls.galaxyDbName = 'VariabilityInfrastructureTestGalaxyDB.db'
        if os.path.exists(cls.galaxyDbName):
            os.unlink(cls.galaxyDbname)

        makeGalTestDB(filename=cls.galaxyDbName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.starDbName):
            os.unlink(cls.starDbName)

        if os.path.exists(cls.galaxyDbName):
            os.unlink(cls.galaxyDbName)

    def setUp(self):
        self.starDB = myTestStars(database=self.starDbName)
        self.galaxyDB = myTestGals(database=self.galaxyDbName)

    def tearDown(self):
        del self.starDB
        del self.galaxyDB


    def testStellarVariabilityInfrastructure(self):
        """
        Test that the variability design was correctly implemented
        in the case of stars
        """

        outputs = ['id', 'test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y']
        baseline = StellarBaselineCatalogClass(self.starDB, column_outputs=outputs)
        variable = StellarVariabilityCatalogClass(self.starDB, column_outputs=outputs)

        for bb, vv in zip(baseline.iter_catalog(), variable.iter_catalog()):
            self.assertEqual(bb[0], vv[0])
            self.assertAlmostEqual(vv[1]-bb[1], 4.0, 10)
            self.assertAlmostEqual(vv[2]-bb[2], 10.0, 10)
            self.assertAlmostEqual(vv[3]-bb[3], 20.0, 10)
            self.assertAlmostEqual(vv[4], bb[4], 10)
            self.assertAlmostEqual(vv[5], bb[5], 10)
            self.assertAlmostEqual(vv[6], bb[6], 10)


    def testGalaxyVariabilityInfrastructure(self):
        """
        Test that the variability design was correctly implemented in
        the case of galaxies
        """

        outputs = ['id',
                   'test_u', 'test_g', 'test_r', 'test_i', 'test_z', 'test_y',
                   'test_bulge_u', 'test_bulge_g', 'test_bulge_r', 'test_bulge_i',
                   'test_bulge_z', 'test_bulge_y',
                   'test_disk_u', 'test_disk_g', 'test_disk_r', 'test_disk_i',
                   'test_disk_z', 'test_disk_y',
                   'test_agn_u', 'test_agn_g', 'test_agn_r',
                   'test_agn_i', 'test_agn_z','test_agn_y']

        baseline = GalaxyBaselineCatalogClass(self.galaxyDB, column_outputs=outputs)
        variable = GalaxyVariabilityCatalogClass(self.galaxyDB, column_outputs=outputs)

        phot = PhotometryGalaxies()

        variable_indices = [19, 20, 21, 7, 16]

        for bb, vv in zip(baseline.iter_catalog(), variable.iter_catalog()):
            self.assertEqual(bb[0], vv[0])

            #test that the variable components are altered
            #the way they ought to be
            self.assertAlmostEqual(bb[19]+1.0, vv[19], 10)
            self.assertAlmostEqual(bb[20]+2.0, vv[20], 10)
            self.assertAlmostEqual(bb[21]+3.0, vv[21], 10)
            self.assertAlmostEqual(bb[7]+4.0, vv[7], 10)
            self.assertAlmostEqual(bb[16]+5.0, vv[16], 10)

            #test that the components which do not vary are equal
            for ix in range(7,25):
                if ix not in variable_indices:
                    self.assertAlmostEqual(bb[ix], vv[ix], 10)

            #test that the total magnitudes are correctly calculated
            for ix in range(6):

                self.assertAlmostEqual(bb[ix+1],
                                       phot.sum_magnitudes(bulge=bb[7+ix],
                                                           disk=bb[13+ix],
                                                           agn=bb[19+ix]),
                                       10)

                self.assertAlmostEqual(vv[ix+1],
                                       phot.sum_magnitudes(bulge=vv[7+ix],
                                                           disk=vv[13+ix],
                                                           agn=vv[19+ix]),
                                       10)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(VariabilityDesignTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)
