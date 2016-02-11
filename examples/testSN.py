"""
SNObject_tests:
A Class containing tests to check crictical functionality for SNObject.py

The following functionality is tested:

    - SED (flambda) for unextincted SEDs in SNCosmo and SNObject
    - SED (flambda) for MW extincted SEDs in SNCosmo and SNObject (independent
        implementations of extinction using OD94 model.)
    - Band Flux for extincted SED in r Band
    - Band Mag for extincted SED in r Band
    - rectification of SED: while the SEDs from the model can go negative
      resulting in negative model fluxes and adus. Test that if rectifySEDs in
      on such situations are avoided
SNIaCatalog_tests:
A Class containing tests to check crictical functionality for SNIaCatalog 
"""
import os
import sqlite3
import numpy as np
import matplotlib.pyplot as plt
import unittest

# Lsst Sims Dependencies
import lsst.utils.tests as utilsTests
from lsst.sims.photUtils import Bandpass
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import ObservationMetaData
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters
from lsst.sims.utils import spatiallySample_obsmetadata as sample_obsmetadata
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
import eups

# Routines Being Tested
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.mixins import SNIaCatalog

# External packages used
# import pandas as pd
from pandas.util.testing import assert_frame_equal
import sncosmo
import astropy


class SNObject_tests(unittest.TestCase):

    def setUp(self):
        """
        Setup tests
        SN_blank: A SNObject with no MW extinction
        """

        from astropy.config import get_config_dir

        mydir = get_config_dir()
        print '==============================='
        print '==============================='
        print (mydir)
        print '==============================='
        print '==============================='
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

        self.lsstBandPass = BandpassDict.loadTotalBandpassesFromFiles()
        self.photParams = PhotometricParameters()
        self.SNCosmoBP = sncosmo.Bandpass(wave=self.lsstBandPass['r'].wavelen,
                                          trans=self.lsstBandPass['r'].sb,
                                          wave_unit=astropy.units.Unit('nm'),
                                          name='lsst_r')

    def tearDown(self):
        pass

    def test_SNstatenotEmpty(self):
        """
        Check that the state of SNObject, stored in self.SNstate has valid
        entries for all keys and does not contain keys with None type Values.
        """
        myDict = self.SN_extincted.SNstate
        for key in myDict.keys():
            assert myDict[key] is not None
                
    def test_rectifiedSED(self):
        """
        Check for an extreme case that the SN seds are being rectified. This is
        done by setting up an extreme case where there will be negative seds, and
        checking that this is indeed the case, and checking that they are not
        negative if rectified.
        """
        
        snobj = SNObject(ra=30., dec=-60., source='salt2')
        snobj.set(z=0.96, t0=self.mjdobs, x1=-3., x0=1.8e-6)
        snobj.rectifySED  = False
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
            print(time-self.mjdobs)
            sed = snobj.SNObjectSED(time=time,
                                    bandpass=self.lsstBandPass['r'])
            assert not sed.calcADU(bandpass=self.lsstBandPass['r'],
                                   photParams=self.photParams) < 0.
            assert not any(sed.flambda < 0.)

    def test_sourceSED(self):
        
        snobject =  SNObject()

        # This is the SNCosmo Source
        SNCosmoSource = sncosmo.get_source(snobject.source.name)

        # Set up a parameter dict
        params = {'x1': -3.0, 'x0': 1.0, 'c':0.}
        SNCosmoSource.set(**params)
        snobject.set(**params)

        rest_wave = np.arange(2000., 9000., 100.)
        for phase in np.arange(-20., 50., 1.0):
            if snobject.rectifySED:
                assert not any(snobject.SNObjectSourceSED(time=phase).flambda < 0.)
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
                                               time=times,  zpsys='ab',
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
                                              time=times,  magsys='ab')
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

class SNIaCatalog_tests(unittest.TestCase):

    @classmethod
    def setUpClass(cls):


        # Set directory where scratch work will be done
        cls.madeScratchDir = False
        cls.scratchDir = 'scratchSpace'

        # Setup a directory in which test data will be made
        if not os.path.exists(cls.scratchDir):
            os.makedirs(cls.scratchDir)
            cls.madeScratchDir = True

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
        cls.GalaxyPositionSamps = sample_obsmetadata(
                obsmetadata=cls.obsMetaDataforCat, size=cls.size)

        # Create a galaxy Table overlapping with the obsMetaData Spatial Bounds
        # using positions from the samples above and a database name given by
        # self.dbname
        vals = cls._createFakeGalaxyDB()
        with open('valsFromTest.dat', 'w') as f:
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
                       ('raJ2000','raJ2000 * PI()/ 180. '),
                       ('decJ2000','decJ2000 * PI()/ 180.'),
                       ('redshift', 'redshift')]

	# class galCopy(InstanceCatalog):
        #   column_outputs = ['id', 'raJ2000', 'decJ2000', 'redshift']
	#   override_formats = {'raJ2000': '%8e', 'decJ2000': '%8e'}

	cls.galDB = MyGalaxyCatalog(database=cls.dbname)
	# cls.galphot = galCopy(db_obj=cls.galDB,
        #  		       obs_metadata=cls.obsMetaDataforCat)
	# cls.galPhotFname = os.path.join(cls.scratchDir, 'gals.dat')
	# cls.galphot.write_catalog(cls.galPhotFname)

        # Generate a set of Observation MetaData Outputs that overlap
        # the galaxies in space
        opsimPath = os.path.join(eups.productDir('sims_data'),'OpSimData')
        opsimDB = os.path.join(opsimPath,'opsimblitz1_1133_sqlite.db')

        generator = ObservationMetaDataGenerator()
        cls.obsMetaDataResults = generator.getObservationMetaData(limit=100,
                                                    fieldRA=(5.0, 8.0), 
                                                    fieldDec=(-85.,-60.),
                                                    expMJD=(49300., 49400.),
                                                    boundLength=0.15,
                                                    boundType='circle')

        # cls.obsMetaDataResults has obsMetaData corresponding to 15 pointings
        # This is tested in test_obsMetaDataGeneration 

#        v = zip(*map(cls.coords, cls.obsMetaDataResults))

#        fig2, ax2 = plt.subplots()
#        ax2.plot(vals[0][:1000], vals[1][: 1000], '.')
#        ax2.plot(v[0], v[1], 'ko', markersize=8)
#        ax2.axhline(-np.pi, color='k', lw=2)
#        ax2.axhline(np.pi, color='k', lw=2)
#        ax2.axvline(0., color='k', lw=2.)
#        ax2.axvline(2. * np.pi, color='k', lw=2.)
#        fig2.savefig(os.path.join(cls.scratchDir, 'matchPointings.pdf'))



        #print 'cls.obsMetaDataforCat'
        #print cls.obsMetaDataforCat.summary

        #print 'obsMetaDataResults'

#         obsMetaDataList = []
#         for obsMetaData in tmpobsMetaDataResults:
#             obsMetaDataList.append(ObservationMetaData(boundType='circle',
#                                    boundLength=np.degrees(0.05),
#                                    unrefractedRA=np.degrees(0.13),
#                                    unrefractedDec=np.degrees(-1.2),
#                                    bandpassName=obsMetaData.bandpass,
#                                    mjd=obsMetaData.mjd))

        # cls.obsMetaDataResults = tmpobsMetaDataResults# pobsMetaDataList

        # self.catalogList = self._writeManySNCatalogs()
        sncatalog = SNIaCatalog(db_obj=cls.galDB,
                                obs_metadata=cls.obsMetaDataResults[12],
                                column_outputs=['t0', 'flux_u', 'flux_g', \
                                                'flux_r', 'flux_i', 'flux_z',\
                                                'flux_y', 'mag_u', 'mag_g',\
                                                'mag_r', 'mag_i', 'mag_z', \
                                                'mag_y', 'adu_u', 'adu_g',\
                                                'adu_r', 'adu_i', 'adu_z', \
                                                'adu_y','mwebv'])
	sncatalog.suppressDimSN = True
        sncatalog.midSurveyTime = sncatalog.mjdobs - 20.
	sncatalog.snFrequency = 1.0
        cls.fullCatalog = cls.scratchDir + '/testSNCatalogTest.dat'
        sncatalog.write_catalog(cls.fullCatalog)

        # Create a SNCatalog based on GalDB, and having times of explosions
        #     overlapping the times in obsMetaData
        cls.fnameList = cls._writeManySNCatalogs(cls.obsMetaDataResults)

    def test_writingfullCatalog(self):
        """
        Check that a full catalog of SN has more than one line
        """

        with open(self.fullCatalog, 'r') as f:
            numLines =  sum(1 for _ in f)

        self.assertGreater(numLines, 1)

    @staticmethod
    def buildLCfromInstanceCatFilenames(fnamelist):
        # External packages used
        import pandas as pd
        from pandas.util.testing import assert_frame_equal
        dfs = []
        _ = map(lambda x: dfs.append(pd.read_csv(x, index_col=None, sep=', ')),
                fnamelist)
        all_lcsDumped = pd.concat(dfs)
        all_lcsDumped.rename(columns={'#snid': 'snid'}, inplace=True)
        all_lcsDumped['snid'] = all_lcsDumped['snid'].astype(int)
        lcs = all_lcsDumped.groupby('snid')

        return lcs

    # Skip the following test using the command below if we integrate into
    # tests before pandas is in
    #@unittest.skip('depends on pandas')
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
        for key in lcs.groups.keys():
            df = lcs.get_group(key)
            for prop in props:
                print(s.format(len(df), df.snid.iloc[0]) + prop)
                np.testing.assert_equal(len(df[prop].unique()), 1)

    # Skip the following test using the command below if we integrate into
    # tests before pandas is in
    #@unittest.skip('depends on  pandas')
    def test_redrawingCatalog(self): 
        """
        test that drawing the same catalog
        """
        from random import shuffle
        import copy

        test_description = 'Compare second draws of catalog to initial draw'
        obsMetaDataResults = copy.deepcopy(self.obsMetaDataResults)
        shuffle(obsMetaDataResults)
        fnameList = self._writeManySNCatalogs(obsMetaDataResults,
                                              suffix='.v2.dat')

        newlcs = self.buildLCfromInstanceCatFilenames(fnameList)
        oldlcs = self.buildLCfromInstanceCatFilenames(self.fnameList)


        for key in oldlcs.groups.keys():
            df_old = oldlcs.get_group(key)
            df_old.sort_values(['time', 'band'], inplace=True)
            df_new = newlcs.get_group(key)
            df_new.sort(['time', 'band'], inplace=True)
            s = "Testing equality for SNID {0:8d} with {1:2d} datapoints" 
            print(s.format(df_new.snid.iloc[0], len(df_old)))
            assert_frame_equal(df_new, df_old)


    @classmethod
    def tearDownClass(cls):
        cls.cleanDB(cls.dbname)
        # for fname in cls.fnameList:
        #   os.remove(fname)
        # os.remove(cls.galPhotFname)
        # os.remove(cls.fullCatalog)
        # If scratch directory was created remove it
        #if cls.madeScratchDir:
        #    os.rmdir(cls.scratchDir)

    def test_obsMetaDataGeneration(self):

        numObs = len(self.obsMetaDataResults)
        self.assertEqual(numObs, 15)



    @staticmethod
    def coords(x): 
            return np.radians(x.summary['unrefractedRA']),\
                np.radians(x.summary['unrefractedDec'])


    @staticmethod
    def cleanDB(dbname, verbose=True):
        '''
        Deletes the database dbname from the disk.
        Parameters
        ----------
        dbname: string, mandatory
            name (abs path) of the database to be deleted
        verbose: Bool, optional, defaults to True
    
        '''
    
        if os.path.exists(dbname):
            if verbose:
                print "deleting database ", dbname
            os.unlink(dbname)
        else:
            if verbose:
                print 'database ', dbname, ' does not exist'

    @classmethod
    def _writeManySNCatalogs(cls, obsMetaDataResults, suffix=''):

        
        fnameList = []
        for obsindex, obsMetaData in enumerate(obsMetaDataResults):

            bandpass =  obsMetaData.bandpass
            cols = ['t0', 'mwebv', 'time', 'band', 'flux', 'flux_err',\
                    'mag', 'mag_err', 'cosmologicalDistanceModulus']
            newCatalog = SNIaCatalog(db_obj=cls.galDB, obs_metadata=obsMetaData,
                                     column_outputs=cols)
            newCatalog.midSurveyTime= 49350
            newCatalog.averageRate = 1.
            newCatalog.suppressDimSN = False
            s = "{0:d}".format(obsindex)
            fname = os.path.join(cls.scratchDir, "SNCatalog_" +  s + suffix)
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
        curs.execute('CREATE TABLE if not exists gals (id INT, raJ2000 FLOAT, decJ2000 FLOAT, redshift FLOAT)')
    
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


#    @staticmethod
#    def sample_obsmetadata(obsmetadata, size=1):
#        '''
#        Sample a square patch on the sphere overlapping obsmetadata
#        field of view by picking the area enclosed in
#        obsmetadata.unrefractedRA \pm obsmetadata.boundLength
#        obsmetadata.unrefractedDec \pm obsmetadata.boundLength
#    
#        Parameters
#        ----------
#        obsmetadata: instance of
#            `sims.catalogs.generation.db.ObservationMetaData`
#    
#        size: integer, optional, defaults to 1
#            number of samples
#    
#    
#        Returns
#        -------
#    
#        tuple of ravals, decvalues
#        '''
#        mydict = obsmetadata.summary
#        phi = np.radians(mydict['pointingRA'])
#        theta = np.radians(mydict['pointingDec'])
#        equalrange = np.radians(mydict['boundLength'])
#        ravals, thetavals = SNIaCatalog_tests.samplePatchOnSphere(phi=phi,
#                                                     theta=theta,
#                                                     delta=equalrange,
#                                                     size=size)
#        return ravals, thetavals
#
#
#
#    @staticmethod
#    def samplePatchOnSphere(phi, theta, delta, size):
#        """
#        Samples of corrdinates in spherical coordinates (\phi,\theta)
#        of length size, uniformly distributed in the region 
#        in 
#        Uniformly distributes samples on a spherical patch between 
#        phi \pm delta and theta \pm delta.
#        
#        Parameters
#        ----------
#        phi: float, mandatory, radians
#            center of the spherical patch in ra with range 
#        theta: float, mandatory, radians
#        delta: float, mandatory, radians
#        size: int, mandatory
#            number of samples
#        """
#        np.random.seed(1)
#        u = np.random.uniform(size=size)
#        v = np.random.uniform(size=size)
#    
#        # phivals = delta * (2. * u - 1) + phi
#        phivals = 2. * delta* u + (phi - delta )
#        phivals = np.where ( phivals >= 0., phivals, phivals + 2. * np.pi)
#        
#        # use conventions in spherical coordinates
#        theta = np.pi/2.0 - theta
#        # thetavals = 2. * delta* v + (theta - delta )
#        # thetavals = np.where ( thetavals < np.pi , thetavals, thetavals - np.pi)
#        # thetavals = np.where ( thetavals > - np.pi , thetavals, thetavals + np.pi)
#        
#        
#        thetamax = theta + delta
#        thetamin = theta - delta
#        # CDF is cos(thetamin) - cos(theta) / cos(thetamin) - cos(thetamax)
#        a = np.cos(thetamin) - np.cos(thetamax)
#        thetavals = np.arccos(-v * a + np.cos(thetamin))
#
#        # Get back to -pi/2 to pi/2 range of decs
#        thetavals = np.pi/2.0 - thetavals 
#
#        return phivals, thetavals

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
                print "deleting database ", dbname
            os.unlink(dbname)
        else:
            if verbose:
                print 'database ', dbname, ' does not exist'

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
        s += "( " + ", ".join(["?"]*len(lst)) + ")"
        return s


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(SNObject_tests)
    suites += unittest.makeSuite(SNIaCatalog_tests)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
