"""
Test the simulation of a SNIa light curve using snObject at specified ra, dec,
and specific SN (SALT2) parameters.
"""
import os
import numpy
import unittest
import matplotlib.pyplot as plt

from lsst.sims.photUtils.Photometry import PhotometryStars, Sed, Bandpass
from lsst.sims.photUtils.Photometry import PhotometryBase
import lsst.utils.tests as utilsTests

from lsst.sims.catUtils.exampleCatalogDefinitions.snObject import SNObject
from lsst.sims.catUtils.exampleCatalogDefinitions.sncat import SNIaCatalog
import utils_for_test as tu

from astropy.units import Unit
from astropy.coordinates import SkyCoord
from sncosmo import Model
import sncosmo

from lsst.sims.catalogs.measures.instance import InstanceCatalog
# from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData

import lsst.sims.catUtils.utils.testUtilsSNe as sq
import sqlite3
# Global Variables for calculating MWEBV through astropy/SNcosmo
dustmaproot = os.getenv('SIMS_DUSTMAPS_DIR')
map_dir = os.path.join(dustmaproot, 'DustMaps')


def _file2lst(fname, i, mjd):
    d = numpy.loadtxt(fname, delimiter=',')
    l = list()
    for i, row in enumerate(d):
        obsid = 'obshist' + str(i)
        lst = [obsid] + [mjd] + row.tolist()
        l.append(lst)
    return l


class testSNIaCatalog(unittest.TestCase):
    """
    Unit tests to test the functionality of the SNIaCatalog class.

    Requires:
    ---------
        connection to LSST database


    Tests:
    ------
        Write out an instance catalog of SNIa
        Find SNIa 'observed' according to obs_metadata associated\
        with an LSST view, catalog in an instance catalog and write to ascii\
        files testData/SNIaCat_i.txt as output.
        testICatOuput:
        testWriteToSQLite:
        testLCFromSQLite:
    """
    mjds = [570123.15 + 3.*i for i in range(2)]

    @classmethod
    def setUpClass(cls):
        # delete previous test db if present
        if os.path.exists('testData/sncat.db'):
            print 'deleting previous database'
            os.unlink('testData/sncat.db')

        mjds = [570123.15 + 3.*i for i in range(2)]
        galDB = CatalogDBObject.from_objid('galaxyTiled')

        for i, myMJD in enumerate(mjds):
            myObsMD = ObservationMetaData(boundType='circle',
                                          boundLength=0.015,
                                          unrefractedRA=5.0,
                                          unrefractedDec=15.0,
                                          bandpassName=
                                          ['u', 'g', 'r', 'i', 'z', 'y'],
                                          mjd=myMJD)
            catalog = SNIaCatalog(db_obj=galDB, obs_metadata=myObsMD)
            fname = "testData/SNIaCat_" + str(i) + ".txt"
            catalog.write_catalog(fname)

    @classmethod
    def tearDownClass(cls):
        pass
        # delete the test db
        # if os.path.exists('testData/sncat.db'):
        #    print 'deleting previous database'
        #    os.unlink('testData/sncat.db')

    def testICatOutput(self):
        """
        Check that the output of the instance catalog SNIaCatalog

        """
        stddata = numpy.loadtxt('testData/SNIaCat_0_std.txt', delimiter=',')
        newdata = numpy.loadtxt('testData/SNIaCat_0.txt', delimiter=',')
        stddata.sort(axis=0)
        newdata.sort(axis=0)
        numpy.testing.assert_allclose(stddata, newdata)

    def testWriteToSQLite(self):
        """
        """
        connection = sqlite3.connect('testData/sncat.db')
        curs = connection.cursor()
        curs.execute('CREATE TABLE if not exists mysncat (id TEXT, mjd FLOAT, snid INT, snra FLOAT, sndec FLOAT, z FLOAT, t0 FLOAT, c FLOAT, x1 FLOAT, x0 FLOAT, mag_u FLOAT, mag_g FLOAT, mag_r FLOAT, mag_i FLOAT, mag_z FLOAT, mag_y FLOAT)')

        for i, myMJD in enumerate(self.mjds):
            fname = "testData/SNIaCat_" + str(i) + ".txt"
            l = _file2lst(fname, i, mjd=myMJD)
            recs = sq.array2dbrecords(l)
            exec_str = sq.insertfromdata(tablename='mysncat', records=recs,
                                         multiple=True)
            curs.executemany(exec_str, recs)
        connection.commit()
        # connection.close()
        # def testLCFromSQLite(self):
        # """
        # """
        # connection = sqlite3.connect('testData/sncat.db')
        # curs = connection.cursor()
        curs.execute('SELECT * FROM mysncat')
        print 'LC In Database: '
        lc = curs.fetchall()
        print type(lc)


class testSNObject(unittest.TestCase):
    """
    Unit tests to test functionality of the SNObject module. The following
    tests are proposed after the initial objects are setup:

    1. testSNObjectnoMWmags_SNCosmo: Compare the lsst band magnitudes computed
    using SNCosmo and the LSST software stack when there is no MW extinction.
    2. testSNObjectMWmags_SNCosmo: Compare the lsst band magnitudes computed
    using SNCosmo and the LSST software stack after the SN Sed has been
    extincted in the MW, according to CCM89 and observed dustmaps.
    3. testSNObject_lc_againstprev: Compare the light curve output by SNObject
    for a certain number of days to an older version of the light curve stored
    on disk in a file 'testData/std_lc.dat'
    """
    def setUp(self):
        """
        Setup the objects used in this test suite class. snObject objs are
        required to have a SN sed namely ra, dec, and SALT model parameters
        x0, x1,c and the redshift of the supernova z

        Store objects:
        self.SNmw : SN with said parameters, exincted in MW
        self.SNnomw : SN with said parameters, unexincted in MW
        self.SNCosmo_nomw :
        self.SNCosmo_mw :
        """
        ra = 204.
        dec = -30.
        # Setup SN object of different kinds we want to check

        # SNObject with ra, dec, whose extinction will be calculated
        # separately

        # Use catsim to calculate SN light curve properties from ra, dec
        # and parameters. Do not use MW, LC properties are then really
        # independent of RA DEC. Store in self.SNmw

        # In order to compare SN light curves between SNCosmo and Catsim
        # we need to set them up separately

        # Object of SNObject class with MW extinction
        SN = SNObject(ra, dec)
        SN.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        self.SNmw = SN

        # SNCosmo Model object with MW extinction. Store in self.SNCosmo_mw
        dust = sncosmo.OD94Dust()
        # dust = sncosmo.CCM89Dust()

        SNCosmo = sncosmo.Model(source='salt2-extended',
                                effects=[dust, dust],
                                effect_names=['host', 'mw'],
                                effect_frames=['rest', 'obs'])

        SNCosmo.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        skycoords = SkyCoord(ra, dec, unit='deg')
        t_mwebv = sncosmo.get_ebv_from_map(skycoords,
                                           mapdir=map_dir,
                                           interpolate=False)
        # Now Set the value of mwebv
        SNCosmo.set(mwebv=t_mwebv)
        self.SNCosmo_mw = SNCosmo

        SNCosmo_nomw = sncosmo.Model(source='salt2-extended',
                                     effects=[dust, dust],
                                     effect_names=['host', 'mw'],
                                     effect_frames=['rest', 'obs'])

        SNCosmo_nomw.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        SNCosmo_nomw.set(mwebv=0.)
        self.SNCosmo_nomw = SNCosmo_nomw

        SNCosmo_float = sncosmo.Model(source='salt2-extended',
                                      effects=[dust, dust],
                                      effect_names=['host', 'mw'],
                                      effect_frames=['rest', 'obs'])

        SNCosmo_float.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        self.SNCosmo_float = SNCosmo_float
        # Setup SNObject with no MW extinction
        SNnoMW = SNObject(ra, dec)
        SNnoMW.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
        SNnoMW.set_MWebv(0.)
        self.SNnoMW = SNnoMW

        # Load LSST sofware bandpass objects for magnitude calculation
        self.bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
        photometry = PhotometryBase()
        photometry.loadBandPassesFromFiles(self.bandPassList)
        self.lsstbands = photometry.bandPassList
        self.times = numpy.arange(-20., 50., 1.0)

        # Load SNCosmo bandpass objects for comparison test
        thisshouldbeNone = tu.getlsstbandpassobjs(self.bandPassList,
                                                  loadcatsim=False, plot=False)
        self.sncosmobands = ['LSST' + band for band in self.bandPassList]

    def testSNObjectnoMWmags_SNCosmo(self):
        """
        magnitude calculation for SNObject without extinction:
        Compare the lsst band magnitudes computed using SNCosmo and the LSST
        software stack. This is done by using the SNObject method bandMags and
        compared with the SNCosmo.Model method bandmags with an SNObject whose
        mwebv attribute is set to zero.

        The stored variables are self.SNnomw and self.SNCosmo_nomw
        """
        lsst = []
        sncosmo = []
        for time in self.times:

            bandMagsfromLSST = self.SNnoMW.bandMags(self.lsstbands, time=time)
            e = [time]
            # e  += bandMagsfromLSST.tolist()
            lsst.append(bandMagsfromLSST.tolist())

            t = time*numpy.ones(len(self.bandPassList))
            t.tolist()
            z = [time]

            self.SNCosmo_nomw.set(mwebv=0.0)
            y = self.SNCosmo_nomw.bandmag(band=self.sncosmobands,
                                          time=t, magsys='ab')
            # z += y.tolist()
            sncosmo.append(y.tolist())
            numpy.testing.assert_allclose(numpy.array(sncosmo),
                                          numpy.array(lsst))

    def testSNmwebv(self):
        """
        The mwebv parameter is found from dustmaps and the coordinates
        (ra, dec) in both the catsim and SNCosmo (astropy) framework. Check
        that these values match to a certain precision
        """
        sncosmoval = self.SNCosmo_mw.parameters[-2]
        lsstval = self.SNmw._mwebv
        numpy.testing.assert_allclose(sncosmoval, lsstval)

    def testSNObjectMWmags_SNCosmowithsamemwebv(self):
        """
        If we use the same mwebv parameter in both the LSST catsim and SNCosmo
        implementation of O'Donnell 94, check for differences in the band
        magnitudes
        """

        lsstval = self.SNmw._mwebv
        self.SNCosmo_float.set(mwebv=lsstval)

        lsst = []
        sncosmo = []
        for time in self.times:

            bandMagsfromLSST = self.SNmw.bandMags(self.lsstbands, time=time)
            e = [time]
            # e  += bandMagsfromLSST.tolist()
            lsst.append(bandMagsfromLSST.tolist())

            t = time*numpy.ones(len(self.bandPassList))
            t.tolist()
            z = [time]
            y = self.SNCosmo_float.bandmag(band=self.sncosmobands, time=t,
                                           magsys='ab')
            # z += y.tolist()
            sncosmo.append(y.tolist())
            numpy.testing.assert_allclose(numpy.array(sncosmo),
                                          numpy.array(lsst))

    def testSNObjectMWmags_SNCosmo(self):
        """
        magnitude calculation for SNObject with extinction:
        Compare the lsst band magnitudes computed using SNCosmo and the LSST
        software stack. This is done by using the SNObject method bandMags and
        compared with the SNCosmo.Model method bandmags with an SNObject whose
        mwebv attribute is calculated in different ways.
        """
        lsst = []
        sncosmo = []
        for time in self.times:

            bandMagsfromLSST = self.SNmw.bandMags(self.lsstbands, time=time)
            e = [time]
            # e  += bandMagsfromLSST.tolist()
            lsst.append(bandMagsfromLSST.tolist())

            t = time*numpy.ones(len(self.bandPassList))
            t.tolist()
            z = [time]
            y = self.SNCosmo_mw.bandmag(band=self.sncosmobands, time=t,
                                        magsys='ab')
            # z += y.tolist()
            sncosmo.append(y.tolist())
            numpy.testing.assert_allclose(numpy.array(sncosmo),
                                          numpy.array(lsst))

    def testSNObject_lc_againstprev(self):
        """
        write out the light curve of a single SN represented as SNObject over a
        number of days, and check that the results are the same as the results
        before stored in testData/std_lc.dat
        """

        l = []
        for time in self.times:
            # for time in numpy.arange(-20., 50., 1.0):
            x = self.SNmw.bandMags(self.lsstbands, time=time)
            # y = SNCosmoSN.bandmag(band=sncosmobands, time=t, magsys='ab')
            e = [time]
            e += x.tolist()
            # e += y.tolist()
            l.append(e)
        header = "time(mjd) u g r i z y"
        lc = numpy.array(l)
        numpy.savetxt('testData/lc.dat', lc, fmt='%10.6f', header=header)

        # The 'testData/standard_lc.dat' file was written using the following
        # for a SN at ra  = 204., -30. ,
        # with x0 = 1.846e-6, x1=0.1, c=0., z=0.2
        # using the SALT2-extended model
        # numpy.savetxt('testData/standard_lc.dat', numpy.array(l),
        #               fmt='%10.6f', header=header)
        std_lc = numpy.loadtxt('testData/standard_lc.dat')
        numpy.testing.assert_allclose(std_lc, lc)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(testSNObject)
    suites += unittest.makeSuite(testSNIaCatalog)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(shouldExit=True)
