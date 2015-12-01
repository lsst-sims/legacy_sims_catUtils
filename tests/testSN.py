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
import os
import sqlite3
import numpy as np
import unittest


# Lsst Sims Dependencies
import lsst.utils.tests as utilsTests
from lsst.sims.photUtils import Bandpass
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
import eups

# Routines Being Tested
from lsst.sims.catUtils.mixins import SNObject
from lsst.sims.catUtils.mixins import SNIaCatalog

# External packages used
import sncosmo
import astropy


class SNObject_tests(unittest.TestCase):

    def setUp(self):
        """
        Setup tests
        SN_blank: A SNObject with no MW extinction


        """

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
        self.SNCosmoBP = sncosmo.Bandpass(wave=self.lsstBandPass['r'].wavelen,
                                          trans=self.lsstBandPass['r'].sb,
                                          wave_unit=astropy.units.Unit('nm'),
                                          name='lsst_r')

    def tearDown(self):
        pass

    def test_ComparebandFluxes2photUtils(self):
        """
        The SNObject.catsimBandFluxes computation uses the sims.photUtils.sed
        band flux computation under the hood. This test makes sure that these
        definitions are in sync
        """

        snobject_r = self.SN_extincted.catsimBandFluxes(
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
        catsim_r = self.SN_extincted.catsimBandFluxes(
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
        catsim_r = self.SN_extincted.catsimBandMags(
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

        print ('-------------------------------')
        print ('Entered SNIaCatalog_tests.Setup')
        print ('-------------------------------')
        # Generate a set of Observation MetaData Outputs that overlap
        #     the galaxies in space
        # Create a SNCatalog based on GalDB, and having times of explosions
        #     overlapping the times in obsMetaData

        # Set directory where scratch work will be done
        cls.madeScratchDir = False
        cls.scratchDir = 'scratchSpace'
        # Setup a directory in which test data will be made
        if not os.path.exists(cls.scratchDir):
            os.makedirs(cls.scratchDir)
            self.madeScratchDir = True

        # ObsMetaData instance with spatial window
        cls.obsMetaDataforCat = ObservationMetaData(boundType='circle',
            boundLength=np.degrees(0.25),
            unrefractedRA=np.degrees(0.13),
            unrefractedDec=np.degrees(-1.2),
            bandpassName=['r'], mjd=49350.)

        # Randomly generate self.size Galaxy positions within the spatial window
        # of obsMetaDataforCat
        cls.dbname = os.path.join(cls.scratchDir, 'galcat.db')
        cls.size = 1000
        cls.GalaxyPositionSamps = cls.sample_obsmetadata(
                obsmetadata=cls.obsMetaDataforCat, size=cls.size)

        # Create a galaxy Table overlapping with the obsMetaData Spatial Bounds
        # using positions from the samples above and a database name given by
        # self.dbname
        vals = cls._createFakeGalaxyDB()

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

	class galCopy(InstanceCatalog):
            column_outputs = ['id', 'raJ2000', 'decJ2000', 'redshift']
	    override_formats = {'raJ2000': '%8e', 'decJ2000': '%8e'}

	cls.galDB = MyGalaxyCatalog(database=cls.dbname)
	cls.galphot = galCopy(db_obj=cls.galDB,
        		       obs_metadata=cls.obsMetaDataforCat)
	fname = os.path.join(cls.scratchDir, 'gals.dat')
	cls.galphot.write_catalog(fname)

        # Setup ObsMetaData Results for SN observation
        opsimPath = os.path.join(eups.productDir('sims_data'),'OpSimData')
        opsimDB = os.path.join(opsimPath,'opsimblitz1_1133_sqlite.db')

        generator = ObservationMetaDataGenerator()
        cls.obsMetaDataResults = generator.getObservationMetaData(limit=100,
                                                    fieldRA=(5.0, 8.0), 
                                                    fieldDec=(-85.,-60.),
                                                    expMJD=(49300., 49400.),
                                                    boundLength=0.015,
                                                    boundType='circle')
        # cls.obsMetaDataResults has obsMetaData corresponding to 15 pointings
        # This is tested in test_obsMetaDataGeneration 


        # self.catalogList = self._writeManySNCatalogs()
        cls.sncatalog = SNIaCatalog(db_obj=cls.galDB,
                                     obs_metadata=cls.obsMetaDataResults[0], 
                                     column_outputs=['t0', 'flux_u', 'flux_g', \
                                                     'flux_r', 'flux_i', 'flux_z',\
                                                     'flux_y', 'mag_u', 'mag_g',\
                                                     'mag_r', 'mag_i', 'mag_z', \
                                                     'mag_y', 'adu_u', 'adu_g',\
                                                     'adu_r', 'adu_i', 'adu_z', \
                                                     'adu_y','mwebv'])
	cls.sncatalog.suppressDimSN = False
        cls.sncatalog.midSurveyTime = cls.sncatalog.mjdobs - 20.
	cls.sncatalog.averageRate = 1.0
        cls.sncatalog.write_catalog(cls.scratchDir + '/testSNCatalog.dat')

    @classmethod
    def tearDown(cls):
        cls.cleanDB(cls.dbname)
        # If scratch directory was created remove it
        if cls.madeScratchDir:
            os.rmdir(scratchDir)

    # def test_obsMetaDataGeneration(self):

    #   numObs = len(self.obsMetaDataResults)
    #    self.assertEqual(numObs, 15)


    def test_GalaxyCatalog(self):

        print "ZHello"
        return 
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
    def _writeManySNCatalogs(cls):

        
        fnameList = []
        for obsindex, obsMetaData in enumerate(cls.obsMetaDataResults):

            print 'iteration number ', obsindex
            # cols = ['t0', 'cosmologicalDistanceModulus', 'mwebv', 'time', \
            #       'band', 'flux', 'flux_err', 'mag', 'mag_err']
            newCatalog = SNIaCatalog(db_obj=cls.galDB, obs_metadata=obsMetaData)
                                     #column_outputs='id')
            newCatalog.midSurveyTime= 49350
            newCatalog.averageRate = 1.
            newCatalog.suppressDimSN = False
            s = "{0:d}".format(obsindex)
            fname = os.path.join(cls.scratchDir, "SNCatalog_" +  s )
            newCatalog.write_catalog(fname)
            fnameList.append(fname)
            print (obsMetaData.mjd)
        return fnameList


    @classmethod
    def _createFakeGalaxyDB(cls, seed=1):
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
    
        np.random.seed(seed)
        #samps = .sample_obsmetadata(ObsMetaData, size=size)
    
        for count in range(size):
            id = 1000000 + count
    
            # Main Database should have values in degrees
            ra = np.degrees(samps[0][count])
            dec = np.degrees(samps[1][count])
            redshift = np.random.uniform()
            row = tuple([id, ra, dec, redshift])
            exec_str = cls.insertfromdata(tablename='gals', records=row,
                                      multiple=False)
            curs.execute(exec_str, row)
    
        conn.commit()
        conn.close()
        return samps


    @staticmethod
    def sample_obsmetadata(obsmetadata, size=1):
        '''
        Sample a square patch on the sphere overlapping obsmetadata
        field of view by picking the area enclosed in
        obsmetadata.unrefractedRA \pm obsmetadata.boundLength
        obsmetadata.unrefractedDec \pm obsmetadata.boundLength
    
        Parameters
        ----------
        obsmetadata: instance of
            `sims.catalogs.generation.db.ObservationMetaData`
    
        size: integer, optional, defaults to 1
            number of samples
    
    
        Returns
        -------
    
        tuple of ravals, decvalues
        '''
        mydict = obsmetadata.summary
        phi = np.radians(mydict['unrefractedRA'])
        theta = np.radians(mydict['unrefractedDec'])
        equalrange = np.radians(mydict['boundLength'])
        ravals, thetavals = SNIaCatalog_tests.samplePatchOnSphere(phi=phi,
                                                     theta=theta,
                                                     delta=equalrange,
                                                     size=size)
        return ravals, thetavals



    @staticmethod
    def samplePatchOnSphere(phi, theta, delta, size):
        """
        Samples of corrdinates in spherical coordinates (\phi,\theta)
        of length size, uniformly distributed in the region 
        in 
        Uniformly distributes samples on a spherical patch between 
        phi \pm delta and theta \pm delta.
        
        Parameters
        ----------
        phi: float, mandatory, radians
            center of the spherical patch in ra with range 
        theta: float, mandatory, radians
        delta: float, mandatory, radians
        size: int, mandatory
            number of samples
        """
        u = np.random.uniform(size=size)
        v = np.random.uniform(size=size)
    
        # phivals = delta * (2. * u - 1) + phi
        phivals = 2. * delta* u + (phi - delta )
        phivals = np.where ( phivals >= 0., phivals, phivals + 2. * np.pi)
        
        # use conventions in spherical coordinates
        theta = np.pi/2.0 - theta
        # thetavals = 2. * delta* v + (theta - delta )
        # thetavals = np.where ( thetavals < np.pi , thetavals, thetavals - np.pi)
        # thetavals = np.where ( thetavals > - np.pi , thetavals, thetavals + np.pi)
        
        
        thetamax = theta + delta
        thetamin = theta - delta
        # CDF is cos(thetamin) - cos(theta) / cos(thetamin) - cos(thetamax)
        a = np.cos(thetamin) - np.cos(thetamax)
        thetavals = np.arccos(-v * a + np.cos(thetamin))

        # Get back to -pi/2 to pi/2 range of decs
        thetavals = np.pi/2.0 - thetavals 

        return phivals, thetavals

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
