from __future__ import with_statement
import os
import unittest
import lsst.utils.tests as utilsTests
import numpy as np
import sqlite3
import json
from lsst.utils import getPackageDir
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.mixins import PhotometryStars, PhotometryGalaxies
from lsst.sims.catUtils.mixins import VariabilityStars, VariabilityGalaxies
from lsst.sims.catUtils.utils import TestVariabilityMixin

def makeMflareTable(size=10, **kwargs):
    """
    Make a test database to serve information to the flare test
    """

    np.random.seed(88)

    #a haphazard sample of mdwarf SEDs
    sedFiles = ['m2.0Full.dat', 'm5.1Full.dat', 'm4.9Full.dat']

    #a haphazard sample of mflare light curves
    lcFiles = ['flare_lc_bin3_4.dat', 'flare_lc_bin1_4.dat', 'flare_lc_bin3_3.dat']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE mFlare
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    for i in xrange(size):
        sedFile = sedFiles[np.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyMflare',
           'pars':{'t0':48000.0, 'lcfilename':lcFiles[np.random.randint(0,len(lcFiles))], 'dt':0.00069444418, 'length': 1825}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO mFlare VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeRRlyTable(size=100, **kwargs):
    """
    Make a test database to serve information to the rrlyrae test
    """

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #a haphazard sample of RRLyrae light curves
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

    np.random.seed(32)
    mjDisplacement = (np.random.sample(size)-50.0)*50.0
    for i in xrange(size):
        sedFile = sedFiles[np.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyRRly',
           'pars':{'tStartMjd':48000.0+mjDisplacement[i], 'filename':lcFiles[np.random.randint(0,len(lcFiles))]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO RRly VALUES (%i, '%s', '%s')''' % (i, paramStr,sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeCepheidTable(size=100, **kwargs):
    """
    Make a test database to serve information to the cepheid test
    """

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #a haphazard sample of cepheid light curves
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

    np.random.seed(32)
    periods = np.random.sample(size)*50.0
    mjDisplacement = (np.random.sample(size)-0.5)*50.0
    for i in xrange(size):
        sedFile = sedFiles[np.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyCepheid',
           'pars':{'period':periods[i], 'lcfile':lcFiles[np.random.randint(0,len(lcFiles))], 't0':48000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO cepheid VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeEbTable(size=100, **kwargs):
    """
    Make a test database to serve information to the Eb test
    """

    #a haphazard sample of eclipsing binary light curves
    lcFiles = ['eb_lc/EB.2294.inp', 'eb_lc/EB.1540.inp', 'eb_lc/EB.2801.inp']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE eb
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    np.random.seed(32)
    periods = np.random.sample(size)*50.0
    mjDisplacement = (np.random.sample(size)-0.5)*50.0
    for i in xrange(size):
        sedFile = 'sed_flat_norm.txt'
        varParam = {'varMethodName':'applyEb',
           'pars':{'period':periods[i], 'lcfile':lcFiles[np.random.randint(0,len(lcFiles))], 't0':48000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO eb VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeMicrolensingTable(size=100, **kwargs):
    """
    Make a test database to serve information to the microlensing test
    """

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #there are two microlensing methods; they should be equivalent
    method = ['applyMicrolensing', 'applyMicrolens']
    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE microlensing
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    np.random.seed(32)
    that = np.random.sample(size)*40.0+40.0
    umin = np.random.sample(size)
    mjDisplacement = np.random.sample(size)*50.0
    for i in xrange(size):
        sedFile = sedFiles[0]
        varParam = {'varMethodName':method[i%len(method)],
           'pars':{'that':that[i], 'umin':umin[i], 't0':52000.0+mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO microlensing VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeBHMicrolensingTable(size=100, **kwargs):
    """
    Make a test database to serve information to the BHmicrolensing test
    """

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #a sample of black hole microlensing light curves that do not repeat time steps
    #(repeating time steps causes the scipy spline interpolation routine to return Nan)
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

    np.random.seed(32)
    mjDisplacement = np.random.sample(size)*5.0*365.25
    for i in xrange(size):
        sedFile = sedFiles[np.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyBHMicrolens',
           'pars':{'filename':lcFiles[np.random.randint(0,len(lcFiles))], 't0':52000.0-mjDisplacement[i]}}
        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO bhmicrolensing VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeAmcvnTable(size=100, **kwargs):
    """
    Make a test database to serve information to the AMCVN test
    """

    #a haphazard sample of white dwarf SEDs
    sedFiles = ['bergeron_He_4750_70.dat_4950', 'bergeron_50000_85.dat_54000']

    conn = sqlite3.connect('VariabilityTestDatabase.db')
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE amcvn
                     (varsimobjid int, variability text, sedfilename text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")

    np.random.seed(32)
    doesBurst = np.random.randint(0,1,size=size)
    burst_freq = np.random.randint(10,150,size=size)
    burst_scale = 115
    amp_burst = np.random.sample(size)*8.0
    color_excess_during_burst = np.random.sample(size)*0.2-0.4
    amplitude = np.random.sample(size)*0.2
    period = np.random.sample(size)*200.0
    mjDisplacement = np.random.sample(size)*50.0
    for i in xrange(size):
        sedFile = sedFiles[np.random.randint(0,len(sedFiles))]
        varParam = {'varMethodName':'applyAmcvn',
           'pars':{'does_burst':int(doesBurst[i]), #have to cast to int from np.int for json
                   'burst_freq':int(burst_freq[i]),
                   'burst_scale':burst_scale,
                   'amp_burst':amp_burst[i],
                   'color_excess_during_burst':color_excess_during_burst[i],
                   'amplitude':amplitude[i],
                   'period':period[i],
                   't0':52000.0+mjDisplacement[i]}}

        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO amcvn VALUES (%i, '%s', '%s')''' % (i, paramStr, sedFile)
        c.execute(qstr)
    conn.commit()
    conn.close()

def makeAgnTable(size=100, **kwargs):
    """
    Make a test database to serve information to the microlensing test
    """

    #a haphazard sample of galaxy SEDs
    sedFiles = ['Exp.31E06.0005Z.spec', 'Inst.79E06.1Z.spec', 'Const.50E07.0005Z.spec']
    method = ['applyAgn']
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

    np.random.seed(32)
    agn_tau = np.random.sample(size)*100.0+100.0
    agn_sfu = np.random.sample(size)*2.0
    agn_sfg = np.random.sample(size)*2.0
    agn_sfr = np.random.sample(size)*2.0
    agn_sfi = np.random.sample(size)*2.0
    agn_sfz = np.random.sample(size)*2.0
    agn_sfy = np.random.sample(size)*2.0
    mjDisplacement = np.random.sample(size)*5.0
    avBulge = np.random.sample(size)*0.5+2.6
    avDisk = np.random.sample(size)*0.5+2.6
    redshift = np.random.sample(size)*0.5
    for i in xrange(size):
        varParam = {'varMethodName':'applyAgn',
           'pars':{'agn_tau':agn_tau[i], 'agn_sfu':agn_sfu[i], 'agn_sfg':agn_sfg[i],
                    'agn_sfr':agn_sfr[i], 'agn_sfi':agn_sfi[i], 'agn_sfz':agn_sfz[i],
                    'agn_sfy':agn_sfy[i], 't0_mjd':48000.0+mjDisplacement[i],
                    'seed':np.random.randint(0,200000)}}

        paramStr = json.dumps(varParam)

        qstr = '''INSERT INTO agn VALUES (%i, %i, %f, %f, %f, '%s', '%s', '%s', '%s')''' % \
               (i, i, avBulge[i], avDisk[i], redshift[i],
               paramStr,
               sedFiles[np.random.randint(0,len(sedFiles))],
               sedFiles[np.random.randint(0,len(sedFiles))],
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

    #a haphazard sample of stellar SEDs
    sedFiles = ['kp10_8750.fits_g35_8950', 'kp03_10500.fits_g45_10600', 'km50_6750.fits_g20_6750']

    #a haphazard sample of cepheid light curves
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

    np.random.seed(32)
    periods = np.random.sample(size)*50.0
    mjDisplacement = (np.random.sample(size)-0.5)*50.0
    for i in xrange(size):
        sedFile = sedFiles[np.random.randint(0,len(sedFiles))]
        if i%3 ==0:
            # just to make sure that Variability mixins no how to andle
            # objects with no variability
            varParam = None
            paramStr = None
        elif i%2 == 0:
            varParam = {'varMethodName':'applyCepheid',
               'pars':{'period':periods[i], 'lcfile':lcFiles[np.random.randint(0,len(lcFiles))], 't0':48000.0+mjDisplacement[i]}}
        else:
            varParam = {'varMethodName':'testVar',
                        'pars':{'period':5.0, 'amplitude':2.0}}

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

class mflareDB(variabilityDB):
    objid = 'mflareTest'
    tableid = 'mFlare'
    objectTypeId = 53

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
    catalog_type = 'stellarVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilename', 'delta_lsst_u']
    default_columns=[('magNorm', 14.0, float)]

class StellarVariabilityCatalogWithTest(InstanceCatalog, PhotometryStars, VariabilityStars, TestVariabilityMixin):
    catalog_type = 'testVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilename', 'delta_lsst_u']
    default_columns=[('magNorm', 14.0, float)]

class OtherVariabilityCatalogWithTest(InstanceCatalog, PhotometryStars, TestVariabilityMixin, VariabilityStars):
    catalog_type = 'otherVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilename', 'delta_lsst_u']
    default_columns=[('magNorm', 14.0, float)]

class GalaxyVariabilityCatalog(InstanceCatalog, PhotometryGalaxies, VariabilityGalaxies):
    catalog_type = 'galaxyVariabilityCatalog'
    column_outputs = ['varsimobjid', 'sedFilenameAgn', 'lsstUdiff', 'delta_uAgn']
    default_columns=[('magNormAgn', 14.0, float),
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
        if os.path.exists('VariabilityTestDatabase.db'):
            os.unlink('VariabilityTestDatabase.db')

    def setUp(self):
        self.scratchSpace = os.path.join(getPackageDir('sims_catUtils'), 'tests')
        self.scratchSpace = os.path.join(self.scratchSpace, 'scratchSpace')
        self.obs_metadata = ObservationMetaData(mjd=52000.0)

    def tearDown(self):
        del self.obs_metadata


    def checkNotEmpty(self, catName):
        """
        Verify that the catalog specified by catName is not empty
        (i.e. that it has more than just a header line)
        """
        with open(catName, 'r') as input_file:
            lines = input_file.readlines()
            self.assertGreater(len(lines), 1)


    def testHybridVariability(self):
        """
        Test that we can generate a catalog which inherits from multiple variability mixins
        (in this case, TestVariability and VariabilityStars).  This is to make sure that
        the register_method and register_class decorators do not mangle inheritance of
        methods from mixins.
        """
        catName = os.path.join(self.scratchSpace, 'var_hybridTestCatalog.dat')
        makeHybridTable()
        myDB = CatalogDBObject.from_objid('hybridTest')
        myCatalog = myDB.getCatalog('testVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

        # make sure order of mixin inheritance does not matter
        catName = os.path.join(self.scratchSpace, 'var_hybridTestCatalog2.dat')
        myCatalog = myDB.getCatalog('otherVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

        # make sure that, if a catalog does not contain a variability method,
        # an error is thrown; verify that it contains the correct error message
        catName = os.path.join(self.scratchSpace, 'var_hybridTestCatalog3.dat')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        with self.assertRaises(RuntimeError) as context:
            myCatalog.write_catalog(catName)

        expectedMessage = "Your InstanceCatalog does not contain a variability method"
        expectedMessage += " corresponding to 'testVar'"
        self.assertTrue(context.exception.message==expectedMessage)

        if os.path.exists(catName):
            os.unlink(catName)


    def testInheritance(self):
        """
        Directly test the contents of the _methodRegistrys for
        StellarVariabilityCatalog and StellarVariabilityCatalogWithTest
        to make sure that method inheritance was handled correctly
        """

        for m1 in StellarVariabilityCatalog._methodRegistry:
            self.assertTrue(m1 in StellarVariabilityCatalogWithTest._methodRegistry)
            self.assertTrue(m1 in OtherVariabilityCatalogWithTest._methodRegistry)

        self.assertTrue('testVar' in StellarVariabilityCatalogWithTest._methodRegistry)
        self.assertTrue('testVar' in OtherVariabilityCatalogWithTest._methodRegistry)
        self.assertFalse('testVar' in StellarVariabilityCatalog._methodRegistry)


    def testMflares(self):
        catName = os.path.join(self.scratchSpace, 'var_mFlareTestCatalog.dat')
        makeMflareTable()
        myDB = CatalogDBObject.from_objid('mflareTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testRRlyrae(self):
        catName = os.path.join(self.scratchSpace, 'var_rrlyTestCatalog.dat')
        makeRRlyTable()
        myDB = CatalogDBObject.from_objid('rrlyTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testCepheids(self):
        catName = os.path.join(self.scratchSpace, 'var_cepheidTestCatalog.dat')
        makeCepheidTable()
        myDB = CatalogDBObject.from_objid('cepheidTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testEb(self):
        catName = os.path.join(self.scratchSpace, 'var_ebTestCatalog.dat')
        makeEbTable()
        myDB = CatalogDBObject.from_objid('ebTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testMicrolensing(self):
        #Note: this test assumes that the parameters for the microlensing variability
        #model occur in a standard varParamStr column in the database.
        #Actually, the current database of microlensing events simply store the variability
        #parameters as independent columns in the database.
        #The varParamStr formalism is how the applyMicrolensing methods are written, however,
        #so that is what we will test.

        catName = os.path.join(self.scratchSpace, 'var_microlensTestCatalog.dat')
        makeMicrolensingTable()
        myDB = CatalogDBObject.from_objid('microlensTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testBHMicrolensing(self):
        #Note: this test assumes that the parameters for the BHmicrolensing variability
        #model occur in a standard varParamStr column in the database.
        #Actually, the current database of BHmicrolensing events simply store the variability
        #parameters as independent columns in the database.
        #The varParamStr formalism is how the applyBHMicrolens method is written, however,
        #so that is what we will test.

        catName = os.path.join(self.scratchSpace, 'var_bhmicrolensTestCatalog.dat')
        makeBHMicrolensingTable()
        myDB = CatalogDBObject.from_objid('bhmicrolensTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testAmcvn(self):
        #Note: this test assumes that the parameters for the Amcvn variability
        #model occur in a standard varParamStr column in the database.
        #Actually, the current database of Amcvn events simply store the variability
        #parameters as independent columns in the database.
        #The varParamStr formalism is how the applyAmcvn method is written, however,
        #so that is what we will test.

        catName = os.path.join(self.scratchSpace, 'var_amcvnTestCatalog.dat')
        makeAmcvnTable()
        myDB = CatalogDBObject.from_objid('amcvnTest')
        myCatalog = myDB.getCatalog('stellarVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

    def testAgn(self):

        catName = os.path.join(self.scratchSpace, 'var_agnTestCatalog.dat')
        makeAgnTable()
        myDB = CatalogDBObject.from_objid('agnTest')
        myCatalog = myDB.getCatalog('galaxyVariabilityCatalog', obs_metadata=self.obs_metadata)
        myCatalog.write_catalog(catName, chunk_size=1000)

        self.checkNotEmpty(catName)

        if os.path.exists(catName):
            os.unlink(catName)

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(VariabilityTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ == "__main__":
    run(True)

