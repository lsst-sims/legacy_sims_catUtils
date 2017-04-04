from __future__ import with_statement
from builtins import range
import os
import unittest
import numpy as np
import sqlite3
import json
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies
from lsst.sims.catUtils.mixins import VariabilityStars, VariabilityGalaxies
from lsst.sims.catUtils.utils import TestVariabilityMixin

from lsst.sims.catUtils.mixins import Variability, reset_agn_lc_cache


def setup_module(module):
    lsst.utils.tests.init()


def makeRRlyTable(size=100, **kwargs):
    """
    Make a test database to serve information to the rrlyrae test
    """

    # a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    # a haphazard sample of RRLyrae light curves
    lcFiles = ['rrly_lc/RRc/959802_per.txt', 'rrly_lc/RRc/1078860_per.txt', 'rrly_lc/RRab/98874_per.txt',
               'rrly_lc/RRab/3879827_per.txt']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE RRly
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    rng = np.random.RandomState(32)
    mjDisplacement = (rng.random_sample(size)-50.0)*50.0
    for i in range(size):
        sedFile = sedFiles[rng.randint(0, len(sedFiles))]
        varParam = {'varMethodName': 'applyRRly',
                    'pars': {'tStartMjd': 48000.0+mjDisplacement[i],
                             'filename': lcFiles[rng.randint(0, len(lcFiles))]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO RRly VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()


def makeCepheidTable(size=100, **kwargs):
    """
    Make a test database to serve information to the cepheid test
    """

    # a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    # a haphazard sample of cepheid light curves
    lcFiles = ['cepheid_lc/classical_longPer_specfile', 'cepheid_lc/classical_medPer_specfile',
               'cepheid_lc/classical_shortPer_specfile', 'cepheid_lc/classical_shortPer_specfile',
               'cepheid_lc/popII_longPer_specfile', 'cepheid_lc/popII_shortPer_specfile']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE cepheid
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    rng = np.random.RandomState(32)
    periods = rng.random_sample(size)*50.0
    mjDisplacement = (rng.random_sample(size)-0.5)*50.0
    for i in range(size):
        sedFile = sedFiles[rng.randint(0, len(sedFiles))]
        varParam = {'varMethodName': 'applyCepheid',
                    'pars': {'period': periods[i], 'lcfile': lcFiles[rng.randint(0, len(lcFiles))],
                             't0': 48000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO cepheid VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()


def makeEbTable(size=100, **kwargs):
    """
    Make a test database to serve information to the Eb test
    """

    # a haphazard sample of eclipsing binary light curves
    lcFiles = ['eb_lc/EB.2294.inp', 'eb_lc/EB.1540.inp', 'eb_lc/EB.2801.inp']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE eb
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    rng = np.random.RandomState(32)
    periods = rng.random_sample(size)*50.0
    mjDisplacement = (rng.random_sample(size)-0.5)*50.0
    for i in range(size):
        sedFile = 'sed_flat_norm.txt'
        varParam = {'varMethodName': 'applyEb',
                    'pars': {'period': periods[i], 'lcfile': lcFiles[rng.randint(0, len(lcFiles))],
                             't0': 48000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO eb VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()


def makeMicrolensingTable(size=100, **kwargs):
    """
    Make a test database to serve information to the microlensing test
    """

    # a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    # there are two microlensing methods; they should be equivalent
    method = ['applyMicrolensing', 'applyMicrolens']
    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE microlensing
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    rng = np.random.RandomState(32)
    that = rng.random_sample(size)*40.0+40.0
    umin = rng.random_sample(size)
    mjDisplacement = rng.random_sample(size)*50.0
    for i in range(size):
        sedFile = sedFiles[0]
        varParam = {'varMethodName': method[i%len(method)],
                    'pars': {'that': that[i], 'umin': umin[i], 't0': 52000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO microlensing VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()


def makeBHMicrolensingTable(size=100, **kwargs):
    """
    Make a test database to serve information to the BHmicrolensing test
    """

    # a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    # a sample of black hole microlensing light curves that do not repeat time steps
    # (repeating time steps causes the scipy spline interpolation routine to return Nan)
    lcFiles = ['microlens/bh_binary_source/lc_14_25_75_8000_0_0.05_316',
               'microlens/bh_binary_source/lc_14_25_4000_8000_0_phi1.09_0.005_100',
               'microlens/bh_binary_source/lc_14_25_75_8000_0_tets2.09_0.005_316']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE bhmicrolensing
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    rng = np.random.RandomState(32)
    mjDisplacement = rng.random_sample(size)*5.0*365.25
    for i in range(size):
        sedFile = sedFiles[rng.randint(0, len(sedFiles))]
        varParam = {'varMethodName': 'applyBHMicrolens',
                    'pars': {'filename': lcFiles[rng.randint(0, len(lcFiles))],
                             't0': 52000.0-mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO bhmicrolensing VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()


def makeAmcvnTable(size=100, **kwargs):
    """
    Make a test database to serve information to the AMCVN test
    """

    # a haphazard sample of white dwarf SEDs
    sedFiles = ['bergeron_He_4750_70.dat_4950', 'bergeron_50000_85.dat_54000']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE amcvn
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    rng = np.random.RandomState(32)
    doesBurst = rng.randint(0, 1, size=size)
    burst_freq = rng.randint(10, 150, size=size)
    burst_scale = 115
    amp_burst = rng.random_sample(size)*8.0
    color_excess_during_burst = rng.random_sample(size)*0.2-0.4
    amplitude = rng.random_sample(size)*0.2
    period = rng.random_sample(size)*200.0
    mjDisplacement = rng.random_sample(size)*50.0
    for i in range(size):
        sedFile = sedFiles[rng.randint(0, len(sedFiles))]
        varParam = {'varMethodName': 'applyAmcvn',
                    'pars': {'does_burst': int(doesBurst[i]),  # have to cast to int from np.int for json
                             'burst_freq': int(burst_freq[i]),
                             'burst_scale': burst_scale,
                             'amp_burst': amp_burst[i],
                             'color_excess_during_burst': color_excess_during_burst[i],
                             'amplitude': amplitude[i],
                             'period': period[i],
                             't0': 52000.0+mjDisplacement[i]}}

        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO amcvn VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()


def makeAgnTable(size=100, **kwargs):
    """
    Make a test database to serve information to the microlensing test
    """

    # a haphazard sample of galaxy SEDs
    sedFiles = ['Exp.31E06.0005Z.spec', 'Inst.79E06.1Z.spec', 'Const.50E07.0005Z.spec']
    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE agn
                     (galid int, varsimobjid int,
                      internalAvBulge real, internalAvDisk real, redshift real,
                      variability text,
                      sedFilenameBulge text, sedFilenameDisk text, sedFilenameAgn text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    rng = np.random.RandomState(32)
    agn_tau = rng.random_sample(size)*100.0+100.0
    agn_sfu = rng.random_sample(size)*2.0
    agn_sfg = rng.random_sample(size)*2.0
    agn_sfr = rng.random_sample(size)*2.0
    agn_sfi = rng.random_sample(size)*2.0
    agn_sfz = rng.random_sample(size)*2.0
    agn_sfy = rng.random_sample(size)*2.0
    mjDisplacement = rng.random_sample(size)*5.0
    avBulge = rng.random_sample(size)*0.5+2.6
    avDisk = rng.random_sample(size)*0.5+2.6
    redshift = rng.random_sample(size)*0.5
    for i in range(size):
        varParam = {'varMethodName': 'applyAgn',
                    'pars': {'agn_tau': agn_tau[i], 'agn_sfu': agn_sfu[i], 'agn_sfg': agn_sfg[i],
                             'agn_sfr': agn_sfr[i], 'agn_sfi': agn_sfi[i], 'agn_sfz': agn_sfz[i],
                             'agn_sfy': agn_sfy[i], 't0_mjd': 48000.0+mjDisplacement[i],
                             'seed': rng.randint(0, 200000)}}

        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO agn VALUES (%i, %i, %f, %f, %f, '%s', '%s', '%s', '%s')''' % \
               (i, i, avBulge[i], avDisk[i], redshift[i],
                paramStr,
                sedFiles[rng.randint(0, len(sedFiles))],
                sedFiles[rng.randint(0, len(sedFiles))],
                'agn.spec')

        c.execute(qstr)
    conn.commit()
    conn.close()


def makeHybridTable(size=100, **kwargs):
    """
    Make a test database that contains a mix of Cepheid variables
    and 'testVar' variables (variables that use the applySineVar
    method defined in the TestVariabilityMixin)
    """

    # a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    # a haphazard sample of cepheid light curves
    lcFiles = ['cepheid_lc/classical_longPer_specfile', 'cepheid_lc/classical_medPer_specfile',
               'cepheid_lc/classical_shortPer_specfile', 'cepheid_lc/classical_shortPer_specfile',
               'cepheid_lc/popII_longPer_specfile', 'cepheid_lc/popII_shortPer_specfile']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE hybrid
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    rng = np.random.RandomState(32)
    periods = rng.random_sample(size)*50.0
    mjDisplacement = (rng.random_sample(size)-0.5)*50.0
    for i in range(size):
        sedFile = sedFiles[rng.randint(0, len(sedFiles))]
        if i%3 == 0:
            # just to make sure that Variability mixins no how to andle
            # objects with no variability
            varParam = None
            paramStr = None
        elif i%2 == 0:
            varParam = {'varMethodName': 'applyCepheid',
                        'pars': {'period': periods[i],
                                 'lcfile': lcFiles[rng.randint(0, len(lcFiles))],
                                 't0': 48000.0+mjDisplacement[i]}}
        else:
            varParam = {'varMethodName': 'testVar',
                        'pars': {'period': rng.random_sample()*100.0, 'amplitude': 2.0}}

        if varParam is not None:
            paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO hybrid VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()


class variabilityDB(CatalogDBObject):
    driver = 'sqlite'
    database = 'VariabilityTestDatabase.db'
    idColKey = 'varsimobjid'
    columns = [('id', 'varsimobjid', int),
               ('sedFilename', 'sedfilename', str, 40),
               ('varParamStr', 'variability', str, 600)]



class hybridDB(variabilityDB):
    objid = 'hybridTest'
    tableid = 'hybrid'
    objectTypeId = 54


class rrlyDB(variabilityDB):
    objid = 'rrlyTest'
    tableid = 'RRly'
    objectTypeId = 55


class cepheidDB(variabilityDB):
    objid = 'cepheidTest'
    tableid = 'cepheid'
    objectTypeId = 56


class ebDB(variabilityDB):
    objid = 'ebTest'
    tableid = 'eb'
    objectTypeId = 57


class microlensDB(variabilityDB):
    objid = 'microlensTest'
    tableid = 'microlensing'
    objectTypeId = 58


class BHmicrolensDB(variabilityDB):
    objid = 'bhmicrolensTest'
    tableid = 'bhmicrolensing'
    objectTypeId = 59


class amcvnDB(variabilityDB):
    objid = 'amcvnTest'
    tableid = 'amcvn'
    objectTypeId = 60


class agnDB(variabilityDB):
    objid = 'agnTest'
    tableid = 'agn'
    objectTypeId = 61


class StellarVariabilityCatalog(InstanceCatalog, PhotometryStars, VariabilityStars):
    catalog_type = __file__ + 'stellarVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilename', 'delta_lsst_u']
    default_columns = [('magNorm', 14.0, float)]


class StellarVariabilityCatalogWithTest(InstanceCatalog, PhotometryStars,
                                        VariabilityStars, TestVariabilityMixin):
    catalog_type = __file__ + 'testVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilename', 'delta_lsst_u', 'varParamStr']
    default_columns = [('magNorm', 14.0, float)]


class OtherVariabilityCatalogWithTest(InstanceCatalog, PhotometryStars,
                                      TestVariabilityMixin, VariabilityStars):
    catalog_type = __file__ + 'other_variability_catalog'
    column_outputs = ['varsimobjid', 'sedFilename', 'delta_lsst_u', 'varParamStr']
    default_columns = [('magNorm', 14.0, float)]


class GalaxyVariabilityCatalog(InstanceCatalog, PhotometryGalaxies, VariabilityGalaxies):
    catalog_type = __file__ + 'galaxyVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilenameAgn', 'lsstUdiff', 'delta_uAgn']
    default_columns = [('magNormAgn', 14.0, float),
                       ('magNormDisk', 14.0, float),
                       ('magNormBulge', 14.0, float)]

    def get_lsstUdiff(self):
        lsstUvar = self.column_by_name('lsst_u')

        bulge = self.column_by_name('uBulge')
        disk = self.column_by_name('uDisk')
        agn = self.column_by_name('uAgn') - self.column_by_name('delta_uAgn')
        lsstU = self.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

        return lsstUvar - lsstU

    def get_agnUdiff(self):
        lsstU = self.column_by_name('uAgn')
        lsstUvar = self.column_by_name('uAgn_var')
        return lsstUvar - lsstU


class VariabilityTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if os.path.exists('VariabilityTestDatabase.db'):
            os.unlink('VariabilityTestDatabase.db')

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists('VariabilityTestDatabase.db'):
            os.unlink('VariabilityTestDatabase.db')

    def setUp(self):
        self.obs_metadata = ObservationMetaData(mjd=52000.0)
        self.scratch_dir = os.path.join(getPackageDir('sims_catUtils'),
                                        'tests', 'scratchSpace')

    def tearDown(self):
        del self.obs_metadata

    def verify_catalogs(self, cat_name):
        """
        Verify that a catalog generated by the unit tests below contains
        the rows it ought to.

        This is done by looking for a corresponding catalog in
        tests/testData and verifying that the two catalogs have identical
        rows.

        Parameters
        ----------
        cat_name is the full path to the test-generated catalog.  The
        comparison catalog will be found in tests/testData
        """

        control_dir = os.path.join(getPackageDir('sims_catUtils'),
                                   'tests', 'testData')
        control_name = cat_name.split('/')[-1]
        control_name = os.path.join(control_dir, control_name)
        with open(control_name, 'r') as control_file:
            control_lines = control_file.readlines()

        with open(cat_name, 'r') as test_file:
            test_lines = test_file.readlines()

        for tt in test_lines:
            self.assertIn(tt, control_lines)

        for cc in control_lines:
            self.assertIn(cc, test_lines)

    def testHybridVariability(self):
        """
        Test that we can generate a catalog which inherits from multiple variability mixins
        (in this case, TestVariability and VariabilityStars).  This is to make sure that
        the register_method and register_class decorators do not mangle inheritance of
        methods from mixins.
        """
        cat_name = os.path.join(self.scratch_dir, 'hybridTestCatalog.dat')
        makeHybridTable()
        myDB = CatalogDBObject.from_objid('hybridTest')
        myCatalog = StellarVariabilityCatalogWithTest(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

        # make sure order of mixin inheritance does not matter
        myCatalog = OtherVariabilityCatalogWithTest(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

        # make sure that, if a catalog does not contain a variability method,
        # an error is thrown; verify that it contains the correct error message
        myCatalog = StellarVariabilityCatalog(myDB, obs_metadata=self.obs_metadata)
        with self.assertRaises(RuntimeError) as context:
            myCatalog.write_catalog(cat_name)

        expectedMessage = "Your InstanceCatalog does not contain a variability method"
        expectedMessage += " corresponding to 'testVar'"
        self.assertEqual(context.exception.args[0], expectedMessage)

        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testInheritance(self):
        """
        Directly test the contents of the _methodRegistrys for
        StellarVariabilityCatalog and StellarVariabilityCatalogWithTest
        to make sure that method inheritance was handled correctly
        """

        for m1 in StellarVariabilityCatalog._methodRegistry:
            self.assertIn(m1, StellarVariabilityCatalogWithTest._methodRegistry)
            self.assertIn(m1, OtherVariabilityCatalogWithTest._methodRegistry)

        self.assertIn('testVar', StellarVariabilityCatalogWithTest._methodRegistry)
        self.assertIn('testVar', OtherVariabilityCatalogWithTest._methodRegistry)
        self.assertNotIn('testVar', StellarVariabilityCatalog._methodRegistry)

    def testRRlyrae(self):
        cat_name = os.path.join(self.scratch_dir, 'rrlyTestCatalog.dat')
        makeRRlyTable()
        myDB = CatalogDBObject.from_objid('rrlyTest')
        myCatalog = StellarVariabilityCatalog(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testCepheids(self):
        cat_name = os.path.join(self.scratch_dir, 'cepheidTestCatalog.dat')
        makeCepheidTable()
        myDB = CatalogDBObject.from_objid('cepheidTest')
        myCatalog = StellarVariabilityCatalog(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testEb(self):
        cat_name = os.path.join(self.scratch_dir, 'ebTestCatalog.dat')
        makeEbTable()
        myDB = CatalogDBObject.from_objid('ebTest')
        myCatalog = StellarVariabilityCatalog(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testMicrolensing(self):
        # Note:  this test assumes that the parameters for the microlensing variability
        # model occur in a standard varParamStr column in the database.
        # Actually, the current database of microlensing events simply store the variability
        # parameters as independent columns in the database.
        # The varParamStr formalism is how the applyMicrolensing methods are written, however,
        # so that is what we will test.

        cat_name = os.path.join(self.scratch_dir, 'microlensTestCatalog.dat')
        makeMicrolensingTable()
        myDB = CatalogDBObject.from_objid('microlensTest')
        myCatalog = StellarVariabilityCatalog(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testBHMicrolensing(self):
        # Note:  this test assumes that the parameters for the BHmicrolensing variability
        # model occur in a standard varParamStr column in the database.
        # Actually, the current database of BHmicrolensing events simply store the variability
        # parameters as independent columns in the database.
        # The varParamStr formalism is how the applyBHMicrolens method is written, however,
        # so that is what we will test.

        cat_name = os.path.join(self.scratch_dir, 'bhmicrolensTestCatalog.dat')
        makeBHMicrolensingTable()
        myDB = CatalogDBObject.from_objid('bhmicrolensTest')
        myCatalog = StellarVariabilityCatalog(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testAmcvn(self):
        # Note:  this test assumes that the parameters for the Amcvn variability
        # model occur in a standard varParamStr column in the database.
        # Actually, the current database of Amcvn events simply store the variability
        # parameters as independent columns in the database.
        # The varParamStr formalism is how the applyAmcvn method is written, however,
        # so that is what we will test.

        cat_name = os.path.join(self.scratch_dir, 'amcvnTestCatalog.dat')
        makeAmcvnTable()
        myDB = CatalogDBObject.from_objid('amcvnTest')
        myCatalog = StellarVariabilityCatalog(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)

    def testAgn(self):

        cat_name = os.path.join(self.scratch_dir, 'agnTestCatalog.dat')
        makeAgnTable()
        myDB = CatalogDBObject.from_objid('agnTest')
        myCatalog = GalaxyVariabilityCatalog(myDB, obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(cat_name, chunk_size=1000)

        self.verify_catalogs(cat_name)
        if os.path.exists(cat_name):
            os.unlink(cat_name)


class AgnCacheTest(unittest.TestCase):

    @classmethod
    def tearDownClass(self):
        sims_clean_up()

    def test_agn_caching(self):
        """
        Test that the light curve caching in applyAgn does not change
        the outcomes of the delta_mag calculations.  We will do this
        by simulating the same AGN twice:  once using caching, once resetting
        the cache after every time step.  We will then verify that the two sets
        of outputs are identical.
        """

        rng = np.random.RandomState(8374)

        var = Variability()

        nn = 5
        seed_list = rng.random_integers(0, 20000, nn)
        toff_list = rng.random_sample(nn)*10000.0+40000.0
        sfz_list = rng.random_sample(nn)*2.0
        sfu_list = rng.random_sample(nn)*2.0
        sfg_list = rng.random_sample(nn)*2.0
        sfr_list = rng.random_sample(nn)*2.0
        sfi_list = rng.random_sample(nn)*2.0
        sfy_list = rng.random_sample(nn)*2.0
        tau_list = rng.random_sample(nn)*20.0+20.0

        param_list = []
        for ix in range(nn):
            params = {'agn_sfu': sfu_list[ix], 'agn_sfg': sfg_list[ix],
                      'agn_sfr': sfr_list[ix], 'agn_sfi': sfi_list[ix],
                      'agn_sfz': sfz_list[ix], 'agn_sfy': sfy_list[ix],
                      't0_mjd': toff_list[ix], 'agn_tau': tau_list[ix],
                      'seed': seed_list[ix]}
            param_list.append(params)

        mjd_list = rng.random_sample(100)*10000.0+50000.0
        mjd_list = np.sort(mjd_list)

        caching_output = []

        for mjd in mjd_list:
            for pp in param_list:
                dd = var.applyAgn(pp, mjd)
                for kk in dd:
                    caching_output.append(dd[kk])

        uncached_output = []
        for mjd in mjd_list:
            for pp in param_list:
                reset_agn_lc_cache()
                dd = var.applyAgn(pp, mjd)
                for kk in dd:
                    uncached_output.append(dd[kk])

        np.testing.assert_array_almost_equal(np.array(caching_output), np.array(uncached_output), decimal=10)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
