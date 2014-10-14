import os
import unittest
import sqlite3
import numpy
import json
from collections import OrderedDict
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData, altAzToRaDec, \
                                             calcObsDefaults, getRotTelPos, Site
from lsst.sims.catUtils.baseCatalogModels import GalaxyBulgeObj, GalaxyDiskObj, GalaxyAgnObj, StarObj
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D, PhoSimCatalogPoint, \
                                                         PhoSimCatalogZPoint
from lsst.sims.coordUtils import AstrometryBase

class SearchReversion(CatalogDBObject):

    def _get_column_query(self, *args, **kwargs):
        return CatalogDBObject._get_column_query(self,*args, **kwargs)

    def _final_pass(self, *args, **kwargs):
        return CatalogDBObject._final_pass(self,*args, **kwargs)

    def query_columns(self, *args, **kwargs):
        return CatalogDBObject.query_columns(self, *args, **kwargs)
                
class testGalaxyBulge(SearchReversion, GalaxyBulgeObj):
    objid = 'phoSimTestBulges'
    objectTypeId = 88
    columns = GalaxyBulgeObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyDisk(SearchReversion, GalaxyDiskObj):
    objid = 'phoSimTestDisks'
    objectTypeId = 89
    columns = GalaxyDiskObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyAgn(SearchReversion, GalaxyAgnObj):
    objid = 'phoSimTestAgn'
    objectTypeId = 90
    columns = GalaxyAgnObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testStars(SearchReversion, StarObj):
    objid = 'phoSimTestStars'
    objectTypeId = 91
    columns = StarObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'galacticAv':
            _to_remove.append(entry)
        elif entry[0] == 'magNorm':
            _to_remove.append(entry)

    for target in _to_remove:
        columns.remove(target)

def makePhoSimTestDB(filename='PhoSimTestDatabase.db', size=1000, seedVal=32, **kwargs):
    """
    Make a test database to serve information to the myTestGals object
    @param size: Number of rows in the database
    @param seedVal: Random seed to use
    """
    
    if os.path.exists(filename):
        os.unlink(filename)
    
    #just an example of some valid SED file names
    galaxy_seds = ['Const.80E07.02z.spec','Inst.80E07.002Z.spec','Burst.19E07.0005Z.spec']
    agn_sed = 'agn.spec'
    star_seds = ['km20_5750.fits_g40_5790','m2.0Full.dat','bergeron_6500_85.dat_6700']
    numpy.random.seed(seedVal)
    
    mjd = 52000.0
    alt = numpy.pi/2.0
    az = 0.0
    band = 'r'
    testSite = Site()
    centerRA, centerDec = altAzToRaDec(alt,az,testSite.longitude,testSite.latitude,mjd)
    rotTel = getRotTelPos(az, centerDec, testSite.latitude, 0.0)

    obsDict = calcObsDefaults(centerRA, centerDec, alt, az, rotTel, mjd, band, 
                 testSite.longitude, testSite.latitude)
    
    obsDict['Opsim_expmjd'] = mjd
    radius = 0.1
    phoSimMetadata = OrderedDict([
                      (k, (obsDict[k],numpy.dtype(type(obsDict[k])))) for k in obsDict])

    obs_metadata = ObservationMetaData(boundType = 'circle', unrefractedRA = numpy.degrees(centerRA),
                                       unrefractedDec = numpy.degrees(centerDec), boundLength = 2.0*radius,
                                       phoSimMetadata=phoSimMetadata)

    conn = sqlite3.connect(filename)
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE galaxy_bulge
                 (galtileid int, galid int, bra real, bdec real, ra real, dec real, magnorm_bulge real,
                 sedname_bulge text, a_b real, b_b real, pa_bulge real, bulge_n int,
                 ext_model_b text, av_b real, rv_b real, u_ab real, g_ab real, r_ab real, i_ab real,
                 z_ab real, y_ab real, redshift real)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating galaxy_bulge table.")

    try:
        c.execute('''CREATE TABLE galaxy
                  (galtileid int, galid int, dra real, ddec real, ra real, dec real,
                  magnorm_disk real, sedname_disk text, a_d real, b_d real, pa_disk real,
                  disk_n int, ext_model_d text, av_d real, rv_d real, u_ab real, 
                  g_ab real, r_ab real, i_ab real, z_ab real, y_ab real, redshift real)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating galaxy table.")

    try:
        c.execute('''CREATE TABLE galaxy_agn
                  (galtileid int, galid int, agnra real, agndec real, ra real, dec real,
                  magnorm_agn real, sedname_agn text, varParamStr text, u_ab real, 
                  g_ab real, r_ab real, i_ab real, z_ab real, y_ab real, redshift real)''')
    except:
        raise RuntimeError("Error creating galaxy_agn table.")

    try:
        c.execute('''CREATE TABLE starsALL_forceseek
                  (simobjid int, ra real, decl real, gal_l real, gal_b real, magNorm real,
                  mudecl real, mura real, galacticAv real, vrad real, varParamStar text, sedFilename text, parallax real)''')
    except:
        raise RuntimeError("Error creating starsALL_forceseek table.")

    rr = numpy.random.sample(size)*numpy.radians(radius)
    theta = numpy.random.sample(size)*2.0*numpy.pi
    ra = numpy.degrees(centerRA + rr*numpy.cos(theta))
    dec = numpy.degrees(centerDec + rr*numpy.sin(theta))
    
    bra = numpy.radians(ra+numpy.random.sample(size)*0.01*radius)
    bdec = numpy.radians(dec+numpy.random.sample(size)*0.01*radius)
    dra = numpy.radians(ra + numpy.random.sample(size)*0.01*radius)
    ddec = numpy.radians(dec + numpy.random.sample(size)*0.01*radius)
    agnra = numpy.radians(ra + numpy.random.sample(size)*0.01*radius)
    agndec = numpy.radians(dec + numpy.random.sample(size)*0.01*radius)
    
    magnorm_bulge = numpy.random.sample(size)*4.0 + 17.0
    magnorm_disk = numpy.random.sample(size)*5.0 + 17.0
    magnorm_agn = numpy.random.sample(size)*5.0 + 17.0
    a_b = numpy.random.sample(size)*0.2
    b_b = numpy.random.sample(size)*0.2
    a_d = numpy.random.sample(size)*0.5
    b_d = numpy.random.sample(size)*0.5

    pa_bulge = numpy.random.sample(size)*360.0
    pa_disk = numpy.random.sample(size)*360.0

    av_b = numpy.random.sample(size)*0.4
    av_d = numpy.random.sample(size)*0.4
    rv_b = numpy.random.sample(size)*0.1 + 3.0
    rv_d = numpy.random.sample(size)*0.1 + 3.0

    u_ab = numpy.random.sample(size)*4.0 + 17.0
    g_ab = numpy.random.sample(size)*4.0 + 17.0
    r_ab = numpy.random.sample(size)*4.0 + 17.0
    i_ab = numpy.random.sample(size)*4.0 + 17.0
    z_ab = numpy.random.sample(size)*4.0 + 17.0
    y_ab = numpy.random.sample(size)*4.0 +17.0
    redshift = numpy.random.sample(size)*2.0
    
    t0_mjd = numpy.random.sample(size)*10.0+mjd
    agn_tau = numpy.random.sample(size)*1000.0 + 1000.0
    agnSeed = numpy.random.random_integers(low=2, high=4000, size=size)
    agn_sfu = numpy.random.sample(size)
    agn_sfg = numpy.random.sample(size)
    agn_sfr = numpy.random.sample(size)
    agn_sfi = numpy.random.sample(size)
    agn_sfz = numpy.random.sample(size)
    agn_sfy = numpy.random.sample(size)
    
    rrStar = numpy.random.sample(size)*numpy.radians(radius)
    thetaStar = numpy.random.sample(size)*2.0*numpy.pi
    raStar = centerRA + rrStar*numpy.cos(thetaStar)
    decStar = centerDec + rrStar*numpy.sin(thetaStar)
    gal_l, gal_b = AstrometryBase.equatorialToGalactic(raStar, decStar)
    
    raStar = numpy.degrees(raStar)
    decStar = numpy.degrees(decStar)
    gal_l = numpy.degrees(gal_l)
    gal_b = numpy.degrees(gal_b)
    
    magnormStar = numpy.random.sample(size)*4.0 + 17.0
    mudecl = numpy.random.sample(size)*0.0001
    mura = numpy.random.sample(size)*0.0001
    galacticAv = numpy.random.sample(size)*0.05*3.1
    vrad = numpy.random.sample(size)*1.0
    parallax = 0.00045+numpy.random.sample(size)*0.00001
    
    print 'ra:',ra
    print 'dec: ',dec
    print 'raStar: ',raStar
    print 'decStar: ',decStar
    
    #just keep filling in

    #in order to really test the phoSim infrastructure, we should make the fake database such
    #that it can be interfaced with the existing catalogmodels in GalaxyModels.py
    
    #okay: what we will do is sub-class GalaxyBulgeObj, GalaxyDiskObj, etc.
    #but then re-direct the query methods to the original CatalogDBObject query
    #methods, rather than the GalaxyTileObj query methods
    
    galtileid = 0
    for i in range(size):

        cmd = '''INSERT INTO galaxy_bulge VALUES (%i, %i, %f, %f, %f, %f, %f,
                     '%s', %f, %f, %f, %i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f)''' %\
                     (i, i, bra[i], bdec[i], ra[i], dec[i], magnorm_bulge[i], galaxy_seds[i%len(galaxy_seds)],
                     a_b[i], b_b[i], pa_bulge[i], 4, 'CCM', av_b[i], rv_b[i], u_ab[i], g_ab[i],
                     r_ab[i], i_ab[i], z_ab[i], y_ab[i], redshift[i])
        c.execute(cmd)

        cmd = '''INSERT INTO galaxy VALUES (%i, %i, %f, %f, %f, %f, %f,
                     '%s', %f, %f, %f, %i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f)''' %\
                     (i, i, dra[i], ddec[i], ra[i], dec[i], magnorm_disk[i],
                     galaxy_seds[i%len(galaxy_seds)], a_d[i], b_d[i], pa_disk[i], 1, 'CCM',
                     av_d[i], rv_d[i], u_ab[i], g_ab[i],
                     r_ab[i], i_ab[i], z_ab[i], y_ab[i], redshift[i])
        c.execute(cmd)

        varParam = {'varMethodName':'applyAgn',
                    'pars':{'agn_tau':agn_tau[i], 't0_mjd':t0_mjd[i],
                    'agn_sfu':agn_sfu[i], 'agn_sfg':agn_sfg[i], 'agn_sfr':agn_sfr[i],
                    'agn_sfi':agn_sfi[i], 'agn_sfz':agn_sfz[i], 'agn_sfy':agn_sfy[i],
                    'seed':int(agnSeed[i])}}

        paramStr = json.dumps(varParam)

        cmd = '''INSERT INTO galaxy_agn VALUES (%i, %i, %f, %f, %f, %f, %f, '%s', '%s',
                                               %f, %f, %f, %f, %f, %f, %f)''' %\
                                               (i, i, agnra[i], agndec[i], ra[i], dec[i],
                                               magnorm_agn[i], agn_sed, paramStr,
                                               u_ab[i], g_ab[i], r_ab[i], i_ab[i],
                                               z_ab[i], y_ab[i], redshift[i])
        c.execute(cmd)


        cmd = '''INSERT INTO starsALL_forceseek VALUES (%i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %s, '%s', %f)''' %\
                  (i, raStar[i], decStar[i], gal_l[i], gal_b[i], magnormStar[i], mudecl[i], mura[i],
                  galacticAv[i], vrad[i], 'NULL', star_seds[i%len(star_seds)], parallax[i])

        c.execute(cmd)

    conn.commit()
    conn.close()
    return obs_metadata


obs_metadata = makePhoSimTestDB(size=10)
bulgeDB = testGalaxyBulge(address='sqlite:///PhoSimTestDatabase.db')
testBulge = PhoSimCatalogSersic2D(bulgeDB, obs_metadata = obs_metadata)
diskDB = testGalaxyDisk(address='sqlite:///PhoSimTestDatabase.db')
testDisk = PhoSimCatalogSersic2D(diskDB, obs_metadata = obs_metadata)
agnDB = testGalaxyAgn(address='sqlite:///PhoSimTestDatabase.db')
testAgn = PhoSimCatalogZPoint(agnDB, obs_metadata = obs_metadata)
starDB = testStars(address='sqlite:///PhoSimTestDatabase.db')
testStar = PhoSimCatalogPoint(starDB, obs_metadata = obs_metadata)

testBulge.write_catalog('phoSim_unittest_example.txt')
testDisk.write_catalog('phoSim_unittest_example.txt', write_header=False, write_mode='a')
testAgn.write_catalog('phoSim_unittest_example.txt', write_header=False, write_mode='a')
testStar.write_catalog('phoSim_unittest_example.txt', write_header=False, write_mode='a')

print agnDB.columns
