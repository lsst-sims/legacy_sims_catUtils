"""
SNObject_tests:
A Class containing tests to check crictical functionality for SNObject.py

The following functionality is tested:

    - SED (flambda) for unextincted SEDs in SNCosmo and SNObject
    - SED (flambda) for MW extincted SEDs in SNCosmo and SNObject (independent
        implementations of extinction using OD94 model.)
    - Band Flux for extincted SED in r Band
    - Band Mag for extincted SED in r Band

SNIaCatalog_tests:
A Class containing tests to check crictical functionality for SNIaCatalog
"""
from __future__ import print_function
from builtins import map
from builtins import str
from builtins import range
import os
import sqlite3
import numpy as np
import unittest

# External packages used
# import pandas as pd
from pandas.util.testing import assert_frame_equal
import sncosmo
import astropy


# Lsst Sims Dependencies
import lsst.utils.tests
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.utils import getPackageDir
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import spatiallySample_obsmetadata as sample_obsmetadata
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catalogs.db import CatalogDBObject, fileDBObject

# Routines Being Tested
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.mixins import SNIaCatalog
from lsst.sims.catUtils.utils import SNIaLightCurveGenerator

# 2016 July 28
# For some reason, the Jenkins slaves used for continuous integration
# cannot properly load the astropy config directories used by sncosmo.
# To prevent this from crashing every build, we will test whether
# the directories can be accessed and, if they cannot, use unittest.skipIf()
# to skip all of the unit tests in this file.
from astropy.config import get_config_dir

_skip_sn_tests = False
try:
    get_config_dir()
except:
    _skip_sn_tests = True


def setup_module(module):
    lsst.utils.tests.init()


@unittest.skipIf(_skip_sn_tests, "cannot properly load astropy config dir")
class SNObject_tests(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()

    def setUp(self):
        """
        Setup tests
        SN_blank: A SNObject with no MW extinction
        """

        mydir = get_config_dir()
        print('===============================')
        print('===============================')
        print (mydir)
        print('===============================')
        print('===============================')
        # A range of wavelengths in Ang
        self.wave = np.arange(3000., 12000., 50.)
        # Equivalent wavelenths in nm
        self.wavenm = self.wave / 10.
        # Time to be used as Peak
        self.mjdobs = 571190

        # Check that we can set up a SED
        # with no extinction
        self.SN_blank = SNObject()
        self.SN_blank.setCoords(ra=30., dec=-60.)
        self.SN_blank.set(z=0.96, t0=571181, x1=2.66, c=0.353, x0=1.796e-6)
        self.SN_blank.set_MWebv(0.)

        self.SN_extincted = SNObject(ra=30., dec=-60.)
        self.SN_extincted.set(z=0.96, t0=571181, x1=2.66, c=0.353,
                              x0=1.796112e-06)

        self.SNCosmoModel = self.SN_extincted.equivalentSNCosmoModel()
        self.rectify_photParams = PhotometricParameters()
        self.lsstBandPass = BandpassDict.loadTotalBandpassesFromFiles()
        self.SNCosmoBP = sncosmo.Bandpass(wave=self.lsstBandPass['r'].wavelen,
                                          trans=self.lsstBandPass['r'].sb,
                                          wave_unit=astropy.units.Unit('nm'),
                                          name='lsst_r')

    def tearDown(self):
        del self.SNCosmoBP
        del self.SN_blank
        del self.SN_extincted

    def test_SNstatenotEmpty(self):
        """
        Check that the state of SNObject, stored in self.SNstate has valid
        entries for all keys and does not contain keys with None type Values.
        """
        myDict = self.SN_extincted.SNstate
        for key in myDict:
            assert myDict[key] is not None

    def test_attributeDefaults(self):
        """
        Check the defaults and the setter properties for rectifySED and
        modelOutSideRange
        """
        snobj = SNObject(ra=30., dec=-60., source='salt2')
        self.assertEqual(snobj.rectifySED, True)
        self.assertEqual(snobj.modelOutSideTemporalRange, 'zero')

        snobj.rectifySED = False
        self.assertFalse(snobj.rectifySED, False)
        self.assertEqual(snobj.modelOutSideTemporalRange, 'zero')

    def test_raisingerror_forunimplementedmodelOutSideRange(self):
        """
        check that correct error is raised if the user tries to assign an
        un-implemented model value to
        `sims.catUtils.supernovae.SNObject.modelOutSideTemporalRange`
        """
        snobj = SNObject(ra=30., dec=-60., source='salt2')
        assert snobj.modelOutSideTemporalRange == 'zero'
        with self.assertRaises(ValueError) as context:
            snobj.modelOutSideTemporalRange = 'False'
        self.assertEqual('Model not implemented, defaulting to zero method\n',
                         context.exception.args[0])

    def test_rectifiedSED(self):
        """
        Check for an extreme case that the SN seds are being rectified. This is
        done by setting up an extreme case where there will be negative seds, and
        checking that this is indeed the case, and checking that they are not
        negative if rectified.
        """

        snobj = SNObject(ra=30., dec=-60., source='salt2')
        snobj.set(z=0.96, t0=self.mjdobs, x1=-3., x0=1.8e-6)
        snobj.rectifySED = False
        times = np.arange(self.mjdobs - 50., self.mjdobs + 150., 1.)
        badTimes = []
        for time in times:
            sed = snobj.SNObjectSED(time=time,
                                    bandpass=self.lsstBandPass['r'])
            if any(sed.flambda < 0.):
                badTimes.append(time)
        # Check that there are negative SEDs
        assert(len(badTimes) > 0)
        snobj.rectifySED = True
        for time in badTimes:
            sed = snobj.SNObjectSED(time=time,
                                    bandpass=self.lsstBandPass['r'])
            self.assertGreaterEqual(sed.calcADU(bandpass=self.lsstBandPass['r'],
                                                photParams=self.rectify_photParams), 0.)
            self.assertFalse(any(sed.flambda < 0.))

    def test_ComparebandFluxes2photUtils(self):
        """
        The SNObject.catsimBandFlux computation uses the sims.photUtils.sed
        band flux computation under the hood. This test makes sure that these
        definitions are in sync
        """

        snobject_r = self.SN_extincted.catsimBandFlux(
            bandpassobject=self.lsstBandPass['r'],
            time=self.mjdobs)

        # `sims.photUtils.Sed`
        sed = self.SN_extincted.SNObjectSED(time=self.mjdobs,
                                            bandpass=self.lsstBandPass['r'])
        sedflux = sed.calcFlux(bandpass=self.lsstBandPass['r'])
        np.testing.assert_allclose(snobject_r, sedflux / 3631.0)

    def test_CompareBandFluxes2SNCosmo(self):
        """
        Compare the r band flux at a particular time computed in SNObject and
        SNCosmo for MW-extincted SEDs. While the underlying sed is obtained
        from SNCosmo the integration with the bandpass is an independent
        calculation in SNCosmo  and catsim
        """

        times = self.mjdobs
        catsim_r = self.SN_extincted.catsimBandFlux(
            bandpassobject=self.lsstBandPass['r'],
            time=times)
        sncosmo_r = self.SNCosmoModel.bandflux(band=self.SNCosmoBP,
                                               time=times, zpsys='ab',
                                               zp=0.)
        np.testing.assert_allclose(sncosmo_r, catsim_r)

    def test_CompareBandMags2SNCosmo(self):
        """
        Compare the r band flux at a particular time computed in SNObject and
        SNCosmo for MW-extincted SEDs. Should work whenever the flux comparison
        above works.
        """
        times = self.mjdobs
        catsim_r = self.SN_extincted.catsimBandMag(
            bandpassobject=self.lsstBandPass['r'],
            time=times)
        sncosmo_r = self.SNCosmoModel.bandmag(band=self.SNCosmoBP,
                                              time=times, magsys='ab')
        np.testing.assert_allclose(sncosmo_r, catsim_r)

    def test_CompareExtinctedSED2SNCosmo(self):
        """
        Compare the extincted SEDS in SNCosmo and SNObject. Slightly more
        non-trivial than comparing unextincted SEDS, as the extinction in
        SNObject uses different code from SNCosmo. However, this is still
        using the same values of MWEBV, rather than reading it off a map.
        """
        SNObjectSED = self.SN_extincted.SNObjectSED(time=self.mjdobs,
                                                    wavelen=self.wavenm)

        SNCosmoSED = self.SNCosmoModel.flux(time=self.mjdobs, wave=self.wave) \
            * 10.

        np.testing.assert_allclose(SNObjectSED.flambda, SNCosmoSED,
                                   rtol=1.0e-7)

    def test_CompareUnextinctedSED2SNCosmo(self):
        """
        Compares the unextincted flux Densities in SNCosmo and SNObject. This
        is mereley a sanity check as SNObject uses SNCosmo under the hood.
        """

        SNCosmoFluxDensity = self.SN_blank.flux(wave=self.wave,
                                                time=self.mjdobs) * 10.

        unextincted_sed = self.SN_blank.SNObjectSED(time=self.mjdobs,
                                                    wavelen=self.wavenm)

        SNObjectFluxDensity = unextincted_sed.flambda
        np.testing.assert_allclose(SNCosmoFluxDensity, SNObjectFluxDensity,
                                   rtol=1.0e-7)

    def test_redshift(self):
        """
        test that the redshift method works as expected by checking that
        if we redshift a SN from its original redshift orig_z to new_z where
        new_z is smaller (larger) than orig_z:
        - 1. x0 increases (decreases)
        - 2. source peak absolute magnitude in BesselB band stays the same
        """
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=70., Om0=0.3)

        orig_z = self.SN_extincted.get('z')
        orig_x0 = self.SN_extincted.get('x0')
        peakabsMag = self.SN_extincted.source_peakabsmag('BessellB', 'AB', cosmo=cosmo)

        lowz = orig_z * 0.5
        highz = orig_z * 2.0

        # Test Case for lower redshift
        self.SN_extincted.redshift(z=lowz, cosmo=cosmo)
        low_x0 = self.SN_extincted.get('x0')
        lowPeakAbsMag = self.SN_extincted.source_peakabsmag('BessellB', 'AB', cosmo=cosmo)

        # Test 1.
        self.assertGreater(low_x0, orig_x0)
        # Test 2.
        self.assertEqual(peakabsMag, lowPeakAbsMag)

        # Test Case for higher redshift
        self.SN_extincted.redshift(z=highz, cosmo=cosmo)
        high_x0 = self.SN_extincted.get('x0')
        HiPeakAbsMag = self.SN_extincted.source_peakabsmag('BessellB', 'AB', cosmo=cosmo)

        # Test 1.
        self.assertLess(high_x0, orig_x0)
        # Test 2.
        self.assertEqual(peakabsMag, HiPeakAbsMag)

    def test_bandFluxErrorWorks(self):
        """
        test that bandflux errors work even if the flux is negative
        """
        times = self.mjdobs

        e = self.SN_extincted.catsimBandFluxError(times,
                                                  self.lsstBandPass['r'],
                                                  m5=24.5, fluxinMaggies=-1.0)
        assert isinstance(e, np.float)
        print(e)
        assert not(np.isinf(e) or np.isnan(e))

    




@unittest.skipIf(_skip_sn_tests, "cannot properly load astropy config dir")
class SNIaCatalog_tests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        # Set directory where scratch work will be done
        cls.scratchDir = os.path.join(getPackageDir('sims_catUtils'), 'tests',
                                      'scratchSpace')

        # ObsMetaData instance with spatial window within which we will
        # put galaxies in a fake galaxy catalog
        cls.obsMetaDataforCat = ObservationMetaData(boundType='circle',
                                                    boundLength=np.degrees(0.25),
                                                    pointingRA=np.degrees(0.13),
                                                    pointingDec=np.degrees(-1.2),
                                                    bandpassName=['r'], mjd=49350.)

        # Randomly generate self.size Galaxy positions within the spatial window
        # of obsMetaDataforCat
        cls.dbname = os.path.join(cls.scratchDir, 'galcat.db')
        cls.size = 1000
        cls.GalaxyPositionSamps = sample_obsmetadata(obsmetadata=cls.obsMetaDataforCat,
                                                     size=cls.size)

        # Create a galaxy Table overlapping with the obsMetaData Spatial Bounds
        # using positions from the samples above and a database name given by
        # self.dbname
        vals = cls._createFakeGalaxyDB()
        cls.valName = os.path.join(cls.scratchDir, 'valsFromTest.dat')
        with open(cls.valName, 'w') as f:
            for i, v in enumerate(vals[0]):
                f.write(str(np.radians(vals[0][i])) + '  ' + str(np.radians(vals[1][i])) + '\n')

        # fig, ax = plt.subplots()
        # ax.plot(vals[0][:1000], vals[1][: 1000], '.')
        # ax.plot([0.13], [-1.2], 'rs', markersize=8)
        # fig.savefig(os.path.join(cls.scratchDir, 'match_galDBPosns.pdf'))

        # Read it into a CatalogDBObject galDB
        class MyGalaxyCatalog(CatalogDBObject):
            '''
            Create a like CatalogDBObject connecting to a local sqlite database
            '''

            objid = 'mytestgals'
            tableid = 'gals'
            idColKey = 'id'
            objectTypeId = 0
            appendint = 10000
            database = cls.dbname
            # dbAddress = './testData/galcat.db'
            raColName = 'raJ2000'
            decColName = 'decJ2000'
            driver = 'sqlite'

            # columns required to convert the ra, dec values in degrees
            # to radians again
            columns = [('id', 'id', int),
                       ('raJ2000', 'raJ2000 * PI()/ 180. '),
                       ('decJ2000', 'decJ2000 * PI()/ 180.'),
                       ('redshift', 'redshift')]

        cls.galDB = MyGalaxyCatalog(database=cls.dbname)

        # Generate a set of Observation MetaData Outputs that overlap
        # the galaxies in space
        opsimPath = os.path.join(getPackageDir('sims_data'), 'OpSimData')
        opsimDB = os.path.join(opsimPath, 'opsimblitz1_1133_sqlite.db')

        generator = ObservationMetaDataGenerator(database=opsimDB)
        cls.obsMetaDataResults = generator.getObservationMetaData(limit=100,
                                                                  fieldRA=(5.0, 8.0),
                                                                  fieldDec=(-85., -60.),
                                                                  expMJD=(49300., 49400.),
                                                                  boundLength=0.15,
                                                                  boundType='circle')

        sncatalog = SNIaCatalog(db_obj=cls.galDB,
                                obs_metadata=cls.obsMetaDataResults[6],
                                column_outputs=['t0', 'flux_u', 'flux_g',
                                                'flux_r', 'flux_i', 'flux_z',
                                                'flux_y', 'mag_u', 'mag_g',
                                                'mag_r', 'mag_i', 'mag_z',
                                                'mag_y', 'adu_u', 'adu_g',
                                                'adu_r', 'adu_i', 'adu_z',
                                                'adu_y', 'mwebv'])
        sncatalog.suppressDimSN = True
        sncatalog.midSurveyTime = sncatalog.mjdobs - 20.
        sncatalog.snFrequency = 1.0
        cls.fullCatalog = cls.scratchDir + '/testSNCatalogTest.dat'
        sncatalog.write_catalog(cls.fullCatalog)

        # Create a SNCatalog based on GalDB, and having times of explosions
        #     overlapping the times in obsMetaData
        cls.fnameList = cls._writeManySNCatalogs(cls.obsMetaDataResults)

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        del cls.galDB
        cls.cleanDB(cls.dbname)
        if os.path.exists(cls.valName):
            os.unlink(cls.valName)

        for fname in cls.fnameList:
            if os.path.exists(fname):
                os.unlink(fname)

        if os.path.exists(cls.fullCatalog):
            os.unlink(cls.fullCatalog)

    def test_writingfullCatalog(self):
        """
        Check that a full catalog of SN has more than one line
        """

        with open(self.fullCatalog, 'r') as f:
            numLines = sum(1 for _ in f)

        self.assertGreater(numLines, 1)

    @staticmethod
    def buildLCfromInstanceCatFilenames(fnamelist):
        # External packages used
        import pandas as pd
        dfs = []
        list(map(lambda x: dfs.append(pd.read_csv(x, index_col=None, sep=', ')),
                 fnamelist))
        all_lcsDumped = pd.concat(dfs)
        all_lcsDumped.rename(columns={'#snid': 'snid'}, inplace=True)
        all_lcsDumped['snid'] = all_lcsDumped['snid'].astype(int)
        lcs = all_lcsDumped.groupby('snid')

        return lcs

    def test_drawReproducibility(self):
        """
        Check that when the same SN (ie. with same snid) is observed with
        different pointings leading to different instance catalogs, the
        values of properties remain the same.
        """
        lcs = self.buildLCfromInstanceCatFilenames(self.fnameList)

        props = ['snid', 'snra', 'sndec', 'z', 'x0', 'x1', 'c',
                 'cosmologicalDistanceModulus', 'mwebv']
        s = "Testing Equality across {0:2d} pointings for reported properties"
        s += " of SN {1:8d} of the property "
        for key in lcs.groups:
            df = lcs.get_group(key)
            for prop in props:
                print(s.format(len(df), df.snid.iloc[0]) + prop)
                np.testing.assert_equal(len(df[prop].unique()), 1)

    def test_redrawingCatalog(self):
        """
        test that drawing the same catalog
        """
        from random import shuffle
        import copy

        obsMetaDataResults = copy.deepcopy(self.obsMetaDataResults)
        shuffle(obsMetaDataResults)
        fnameList = self._writeManySNCatalogs(obsMetaDataResults,
                                              suffix='.v2.dat')

        newlcs = self.buildLCfromInstanceCatFilenames(fnameList)
        oldlcs = self.buildLCfromInstanceCatFilenames(self.fnameList)

        for key in oldlcs.groups:
            df_old = oldlcs.get_group(key)
            df_old.sort_values(['time', 'band'], inplace=True)
            df_new = newlcs.get_group(key)
            df_new.sort_values(['time', 'band'], inplace=True)
            s = "Testing equality for SNID {0:8d} with {1:2d} datapoints"
            print(s.format(df_new.snid.iloc[0], len(df_old)))
            assert_frame_equal(df_new, df_old)

        for fname in fnameList:
            if os.path.exists(fname):
                os.unlink(fname)

    def test_obsMetaDataGeneration(self):

        numObs = len(self.obsMetaDataResults)
        self.assertEqual(numObs, 15)

    @staticmethod
    def coords(x):
        return np.radians(x.summary['unrefractedRA']),\
            np.radians(x.summary['unrefractedDec'])

    @classmethod
    def _writeManySNCatalogs(cls, obsMetaDataResults, suffix=''):

        fnameList = []
        for obsindex, obsMetaData in enumerate(obsMetaDataResults):

            cols = ['t0', 'mwebv', 'time', 'band', 'flux', 'flux_err',
                    'mag', 'mag_err', 'cosmologicalDistanceModulus']
            newCatalog = SNIaCatalog(db_obj=cls.galDB, obs_metadata=obsMetaData,
                                     column_outputs=cols)
            newCatalog.midSurveyTime = 49350
            newCatalog.averageRate = 1.
            newCatalog.suppressDimSN = False
            s = "{0:d}".format(obsindex)
            fname = os.path.join(cls.scratchDir, "SNCatalog_" + s + suffix)
            newCatalog.write_catalog(fname)
            fnameList.append(fname)
        return fnameList

    @classmethod
    def _createFakeGalaxyDB(cls):
        '''
        Create a local sqlite galaxy database having filename dbname with

        variables id, raJ2000, decJ2000 and redshift, having number of
        rows =size, and having overlap with ObsMetaData.

        Parameters
        ----------

        '''
        dbname = cls.dbname
        samps = cls.GalaxyPositionSamps
        size = cls.size
        cls.cleanDB(dbname)
        conn = sqlite3.connect(dbname)
        curs = conn.cursor()
        curs.execute('CREATE TABLE if not exists gals '
                     '(id INT, raJ2000 FLOAT, decJ2000 FLOAT, redshift FLOAT)')

        seed = 1
        np.random.seed(seed)

        for count in range(size):
            id = 1000000 + count

            # Main Database should have values in degrees
            ra = samps[0][count]
            dec = samps[1][count]
            redshift = np.random.uniform()
            row = tuple([id, ra, dec, redshift])
            exec_str = cls.insertfromdata(tablename='gals', records=row,
                                          multiple=False)
            curs.execute(exec_str, row)

        conn.commit()
        conn.close()
        return samps

    @staticmethod
    def cleanDB(dbname, verbose=True):
        """
        Deletes the database dbname from the disk.
        Parameters
        ----------
        dbname: string, mandatory
            name (abs path) of the database to be deleted
        verbose: Bool, optional, defaults to True
        """

        if os.path.exists(dbname):
            if verbose:
                print("deleting database ", dbname)
            os.unlink(dbname)
        else:
            if verbose:
                print('database ', dbname, ' does not exist')

    @staticmethod
    def insertfromdata(tablename, records, multiple=True):
        """
        construct string to insert multiple records into sqlite3 database
        args:
            tablename: str, mandatory
                Name of table in the database.
            records: set of records
            multiple:
        returns:
        """
        if multiple:
            lst = records[0]
        else:
            lst = records
        s = 'INSERT INTO ' + str(tablename) + ' VALUES '
        s += "( " + ", ".join(["?"] * len(lst)) + ")"
        return s


class SNIaLightCurveControlCatalog(SNIaCatalog):
    catalog_type = __file__ + 'sn_ia_lc_cat'
    column_outputs = ['uniqueId', 'flux', 'flux_err', 'redshift']
    _midSurveyTime = 49000.0
    _snFrequency = 0.001


@unittest.skipIf(_skip_sn_tests, "cannot properly load astropy config dir")
class SNIaLightCurveTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        rng = np.random.RandomState(99)
        n_sne = 100
        ra_list = rng.random_sample(n_sne) * 7.0 + 78.0
        dec_list = rng.random_sample(n_sne) * 4.0 - 69.0
        zz_list = rng.random_sample(n_sne) * 1.0 + 0.05

        cls.input_cat_name = os.path.join(getPackageDir("sims_catUtils"), "tests")
        cls.input_cat_name = os.path.join(cls.input_cat_name, "scratchSpace", "sne_input_cat.txt")

        with open(cls.input_cat_name, "w") as output_file:
            for ix in range(n_sne):
                output_file.write("%d;%.12f;%.12f;%.12f;%.12f;%.12f\n"
                                  % (ix + 1, ra_list[ix], dec_list[ix],
                                     np.radians(ra_list[ix]), np.radians(dec_list[ix]),
                                     zz_list[ix]))

        dtype = np.dtype([('id', np.int),
                          ('raDeg', np.float), ('decDeg', np.float),
                          ('raJ2000', np.float), ('decJ2000', np.float),
                          ('redshift', np.float)])

        cls.db = fileDBObject(cls.input_cat_name, delimiter=';',
                              runtable='test', dtype=dtype,
                              idColKey='id')

        cls.db.raColName = 'raDeg'
        cls.db.decColName = 'decDeg'
        cls.db.objectTypeId = 873

        cls.opsimDb = os.path.join(getPackageDir("sims_data"), "OpSimData")
        cls.opsimDb = os.path.join(cls.opsimDb, "opsimblitz1_1133_sqlite.db")

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.input_cat_name):
            os.unlink(cls.input_cat_name)

    def test_sne_light_curves(self):
        """
        Generate some super nova light curves.  Verify that they come up with the same
        magnitudes and uncertainties as supernova catalogs.
        """

        gen = SNIaLightCurveGenerator(self.db, self.opsimDb)

        raRange = (78.0, 85.0)
        decRange = (-69.0, -65.0)
        bandpass = 'r'

        pointings = gen.get_pointings(raRange, decRange, bandpass=bandpass)
        gen.sn_universe._midSurveyTime = 49000.0
        gen.sn_universe._snFrequency = 0.001
        self.assertGreater(len(pointings), 1)
        lc_dict, truth = gen.light_curves_from_pointings(pointings)
        self.assertGreater(len(lc_dict), 0)

        for group in pointings:
            self.assertGreater(len(group), 1)
            for obs in group:
                cat = SNIaLightCurveControlCatalog(self.db, obs_metadata=obs)
                for sn in cat.iter_catalog():
                    if sn[1] > 0.0:
                        lc = lc_dict[sn[0]][bandpass]
                        dex = np.argmin(np.abs(lc['mjd'] - obs.mjd.TAI))
                        self.assertLess(np.abs(lc['mjd'][dex] - obs.mjd.TAI), 1.0e-7)
                        self.assertLess(np.abs(lc['flux'][dex] - sn[1]), 1.0e-7)
                        self.assertLess(np.abs(lc['error'][dex] - sn[2]), 1.0e-7)

    def test_sne_light_curves_z_cut(self):
        """
        Generate some super nova light curves.  Add a cutoff in redshift.
        Verify that they come up with the same magnitudes and uncertainties
        as supernova catalogs and that objects with z>z_cutoff are not returned.
        """
        z_cut = 0.9

        gen = SNIaLightCurveGenerator(self.db, self.opsimDb)
        gen.z_cutoff = z_cut

        raRange = (78.0, 85.0)
        decRange = (-69.0, -65.0)
        bandpass = 'r'

        pointings = gen.get_pointings(raRange, decRange, bandpass=bandpass)
        gen.sn_universe._midSurveyTime = 49000.0
        gen.sn_universe._snFrequency = 0.001
        self.assertGreater(len(pointings), 1)
        lc_dict, truth = gen.light_curves_from_pointings(pointings)
        self.assertGreater(len(lc_dict), 0)

        over_z = 0

        for group in pointings:
            self.assertGreater(len(group), 1)
            for obs in group:
                cat = SNIaLightCurveControlCatalog(self.db, obs_metadata=obs)
                for sn in cat.iter_catalog():
                    if sn[1] > 0.0:
                        if sn[3] > z_cut:
                            self.assertNotIn(sn[0], lc_dict)
                            over_z += 1
                        else:
                            lc = lc_dict[sn[0]][bandpass]
                            dex = np.argmin(np.abs(lc['mjd'] - obs.mjd.TAI))
                            self.assertLess(np.abs(lc['mjd'][dex] - obs.mjd.TAI), 1.0e-7)
                            self.assertLess(np.abs(lc['flux'][dex] - sn[1]), 1.0e-7)
                            self.assertLess(np.abs(lc['error'][dex] - sn[2]),
                                            1.0e-7, msg='%e vs %e' % (lc['error'][dex], sn[2]))

        self.assertGreater(over_z, 0)

    def test_sne_multiband_light_curves(self):
        """
        Generate some super nova light curves.  Verify that they come up with the same
        magnitudes and uncertainties as supernova catalogs.  Use multiband light curves.
        """

        gen = SNIaLightCurveGenerator(self.db, self.opsimDb)

        raRange = (78.0, 85.0)
        decRange = (-69.0, -65.0)

        pointings = gen.get_pointings(raRange, decRange, bandpass=('r', 'z'))
        gen.sn_universe._midSurveyTime = 49000.0
        gen.sn_universe._snFrequency = 0.001
        self.assertGreater(len(pointings), 1)
        lc_dict, truth = gen.light_curves_from_pointings(pointings)
        self.assertGreater(len(lc_dict), 0)

        obs_gen = ObservationMetaDataGenerator(database=self.opsimDb, driver='sqlite')
        control_obs_r = obs_gen.getObservationMetaData(fieldRA=raRange, fieldDec=decRange,
                                                       telescopeFilter='r', boundLength=1.75)

        control_obs_z = obs_gen.getObservationMetaData(fieldRA=raRange, fieldDec=decRange,
                                                       telescopeFilter='z', boundLength=1.75)

        self.assertGreater(len(control_obs_r), 0)
        self.assertGreater(len(control_obs_z), 0)

        ct_r = 0
        for obs in control_obs_r:
            cat = SNIaLightCurveControlCatalog(self.db, obs_metadata=obs)
            for sn in cat.iter_catalog():
                if sn[1] > 0.0:
                    ct_r += 1
                    lc = lc_dict[sn[0]]['r']
                    dex = np.argmin(np.abs(lc['mjd'] - obs.mjd.TAI))
                    self.assertLess(np.abs(lc['mjd'][dex] - obs.mjd.TAI), 1.0e-7)
                    self.assertLess(np.abs(lc['flux'][dex] - sn[1]), 1.0e-7)
                    self.assertLess(np.abs(lc['error'][dex] - sn[2]), 1.0e-7)

        self.assertGreater(ct_r, 0)

        ct_z = 0
        for obs in control_obs_z:
            cat = SNIaLightCurveControlCatalog(self.db, obs_metadata=obs)
            for sn in cat.iter_catalog():
                if sn[1] > 0.0:
                    ct_z += 1
                    lc = lc_dict[sn[0]]['z']
                    dex = np.argmin(np.abs(lc['mjd'] - obs.mjd.TAI))
                    self.assertLess(np.abs(lc['mjd'][dex] - obs.mjd.TAI), 1.0e-7)
                    self.assertLess(np.abs(lc['flux'][dex] - sn[1]), 1.0e-7)
                    self.assertLess(np.abs(lc['error'][dex] - sn[2]), 1.0e-7)

        self.assertGreater(ct_z, 0)

    def test_limit_sne_light_curves(self):
        """
        Test that we can limit the number of light curves returned per field of view
        """
        lc_limit = 2
        gen = SNIaLightCurveGenerator(self.db, self.opsimDb)
        gen.sn_universe._midSurveyTime = 49000.0
        gen.sn_universe._snFrequency = 0.001

        raRange = (78.0, 85.0)
        decRange = (-69.0, -65.0)

        pointings = gen.get_pointings(raRange, decRange, bandpass=('r', 'z'))

        control_lc, truth = gen.light_curves_from_pointings(pointings)
        test_lc, truth = gen.light_curves_from_pointings(pointings, lc_per_field=lc_limit)
        self.assertGreater(len(control_lc), len(test_lc))
        self.assertLessEqual(len(test_lc), lc_limit*len(pointings))


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
