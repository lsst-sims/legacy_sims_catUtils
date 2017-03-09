from builtins import range
import numpy as np
import os
import sqlite3
import json
from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters
from lsst.sims.utils import (_raDecFromAltAz, _getRotSkyPos, _getRotTelPos, Site,
                             raDecFromAltAz, haversine, ObservationMetaData)


__all__ = ["calcADUwrapper", "makePhoSimTestDB"]


def calcADUwrapper(sedName=None, magNorm=None, redshift=None, internalAv=None, internalRv=None,
                   galacticAv=None, galacticRv=None, bandpass=None):
    """
    Read in an SED and calculat the number of ADU produced by that SED in a specified bandpass

    Parameters
    ----------
    sedName is a string specifying the file name of the SED

    magNorm is the normalizing magnitude of the SED in the imsimBandpass

    redshift is the redshift of the SED

    internalAv is the Av due to internal dust of the source (if a galaxy)

    internalRv is the Rv due to internal dust of the source (if a galaxy)

    galacticAv is the Av due to Milky Way dust between observer and source

    galacticRv is the Rv due to Milky Way dust between observer and source

    bandpass is an intantiation of Bandpass representing the band in which the ADUs are measured

    Returns
    -------
    A float representing the number of ADUs measured in the bandpass
    """

    imsimband = Bandpass()
    imsimband.imsimBandpass()
    sed = Sed()
    sed.readSED_flambda(sedName)
    fNorm = sed.calcFluxNorm(magNorm, imsimband)
    sed.multiplyFluxNorm(fNorm)
    if internalAv is not None and internalRv is not None:
        if internalAv != 0.0 and internalRv != 0.0:
            a_int, b_int = sed.setupCCMab()
            sed.addCCMDust(a_int, b_int, A_v=internalAv, R_v=internalRv)

    if redshift is not None and redshift != 0.0:
        sed.redshiftSED(redshift, dimming=True)

    a_int, b_int = sed.setupCCMab()
    sed.addCCMDust(a_int, b_int, A_v=galacticAv, R_v=galacticRv)

    adu = sed.calcADU(bandpass, photParams=PhotometricParameters())

    return adu


def makePhoSimTestDB(filename='PhoSimTestDatabase.db', size=1000, seedVal=32, radius=0.1,
                     deltaRA=None, deltaDec=None,
                     bandpass='r', m5=None, seeing=None, **kwargs):
    """
    Make a test database to storing cartoon information for the test phoSim input
    catalog to use.

    The method will return an ObservationMetaData object guaranteed to encompass the
    objects in this database.

    @param [in] filename is a string indicating the name of the DB file to be created

    @param [in] size is the number of objects int he database

    @param [in] seedVal is the seed passed to the random number generator

    @param [in] radius is the radius (in degrees) of the field of view to be returned

    @param [in] bandpass is the bandpas(es) of the observation to be passed to
    ObservationMetaData (optional)

    @param [in] m5 is the m5 value(s) to be passed to ObservationMetaData
    (optional)

    @param [in] seeing is the seeing value(s) in arcseconds to be passed to
    ObservationMetaData (optional)

    @param [in] deltaRA/Dec are numpy arrays that indicate where (in relation to the center
    of the field of view) objects should be placed.  These coordinates are in degrees.  Specifying
    either of these paramters will overwrite size.  If you only specify one of these parameters, the other
    will be set randomly.  These parameters are optional.
    """

    if os.path.exists(filename):
        os.unlink(filename)

    # just an example of some valid SED file names
    galaxy_seds = ['Const.80E07.02Z.spec', 'Inst.80E07.002Z.spec', 'Burst.19E07.0005Z.spec']
    agn_sed = 'agn.spec'
    star_seds = ['km20_5750.fits_g40_5790', 'm2.0Full.dat', 'bergeron_6500_85.dat_6700']

    rng = np.random.RandomState(seedVal)

    if deltaRA is not None and deltaDec is not None:
        if len(deltaRA) != len(deltaDec):
            raise RuntimeError("WARNING in makePhoSimTestDB deltaRA and "
                               "deltaDec have different lengths")

    if deltaRA is not None:
        size = len(deltaRA)
    elif deltaDec is not None:
        size = len(deltaDec)

    # create the ObservationMetaData object
    mjd = 52000.0
    alt = np.pi/2.0
    az = 0.0

    testSite = Site(name='LSST')
    obsTemp = ObservationMetaData(mjd=mjd, site=testSite)
    centerRA, centerDec = _raDecFromAltAz(alt, az, obsTemp)
    rotTel = _getRotTelPos(centerRA, centerDec, obsTemp, 0.0)
    rotSkyPos = _getRotSkyPos(centerRA, centerDec, obsTemp, rotTel)

    obs_metadata = ObservationMetaData(pointingRA=np.degrees(centerRA),
                                       pointingDec=np.degrees(centerDec),
                                       rotSkyPos=np.degrees(rotSkyPos),
                                       bandpassName=bandpass,
                                       mjd=mjd,
                                       boundType = 'circle', boundLength = 2.0*radius,
                                       site=testSite,
                                       m5=m5, seeing=seeing)

    moon_alt = -90.0
    sun_alt = -90.0

    moon_ra, moon_dec = raDecFromAltAz(moon_alt, 0.0, obs_metadata)
    dist2moon = haversine(np.radians(moon_ra), np.radians(moon_dec),
                          obs_metadata._pointingRA, obs_metadata._pointingDec)

    obs_metadata.OpsimMetaData = {'moonra': moon_ra,
                                  'moondec': moon_dec,
                                  'moonalt': moon_alt,
                                  'sunalt': sun_alt,
                                  'dist2moon': dist2moon,
                                  'rottelpos': np.degrees(rotTel)}

    # Now begin building the database.
    # First create the tables.
    conn = sqlite3.connect(filename)
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE galaxy_bulge
                 (galtileid int, galid int, bra real, bdec real, ra real, dec real, magnorm_bulge real,
                 sedname_bulge text, a_b real, b_b real, pa_bulge real, bulge_n int,
                 ext_model_b text, av_b real, rv_b real, u_ab real, g_ab real, r_ab real, i_ab real,
                 z_ab real, y_ab real, redshift real, BulgeHalfLightRadius real)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating galaxy_bulge table.")

    try:
        c.execute('''CREATE TABLE galaxy
                     (galtileid int, galid int, ra real, dec real,
                      bra real, bdec real, dra real, ddec real,
                      agnra real, agndec real,
                      magnorm_bulge, magnorm_disk, magnorm_agn,
                      sedname_bulge text, sedname_disk text, sedname_agn text,
                      varParamStr text,
                      a_b real, b_b real, pa_bulge real, bulge_n int,
                      a_d real, b_d real, pa_disk real, disk_n int,
                      ext_model_b text, av_b real, rv_b real,
                      ext_model_d text, av_d real, rv_d real,
                      u_ab real, g_ab real, r_ab real, i_ab real,
                      z_ab real, y_ab real,
                      redshift real, BulgeHalfLightRadius real, DiskHalfLightRadius real)''')

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
        c.execute('''CREATE TABLE StarAllForceseek
                  (simobjid int, ra real, decl real, magNorm real,
                  mudecl real, mura real, galacticAv real, vrad real, varParamStr text,
                  sedFilename text, parallax real)''')
    except:
        raise RuntimeError("Error creating StarAllForceseek table.")

    # Now generate the data to be stored in the tables.

    rr = rng.random_sample(size)*np.radians(radius)
    theta = rng.random_sample(size)*2.0*np.pi

    if deltaRA is None:
        ra = np.degrees(centerRA + rr*np.cos(theta))
    else:
        ra = np.degrees(centerRA) + deltaRA

    if deltaDec is None:
        dec = np.degrees(centerDec + rr*np.sin(theta))
    else:
        dec = np.degrees(centerDec) + deltaDec

    bra = np.radians(ra+rng.random_sample(size)*0.01*radius)
    bdec = np.radians(dec+rng.random_sample(size)*0.01*radius)
    dra = np.radians(ra + rng.random_sample(size)*0.01*radius)
    ddec = np.radians(dec + rng.random_sample(size)*0.01*radius)
    agnra = np.radians(ra + rng.random_sample(size)*0.01*radius)
    agndec = np.radians(dec + rng.random_sample(size)*0.01*radius)

    magnorm_bulge = rng.random_sample(size)*4.0 + 17.0
    magnorm_disk = rng.random_sample(size)*5.0 + 17.0
    magnorm_agn = rng.random_sample(size)*5.0 + 17.0
    b_b = rng.random_sample(size)*0.2
    a_b = b_b+rng.random_sample(size)*0.05
    b_d = rng.random_sample(size)*0.5
    a_d = b_d+rng.random_sample(size)*0.1

    BulgeHalfLightRadius = rng.random_sample(size)*0.2
    DiskHalfLightRadius = rng.random_sample(size)*0.5

    pa_bulge = rng.random_sample(size)*360.0
    pa_disk = rng.random_sample(size)*360.0

    av_b = rng.random_sample(size)*0.3
    av_d = rng.random_sample(size)*0.3
    rv_b = rng.random_sample(size)*0.1 + 2.0
    rv_d = rng.random_sample(size)*0.1 + 2.0

    u_ab = rng.random_sample(size)*4.0 + 17.0
    g_ab = rng.random_sample(size)*4.0 + 17.0
    r_ab = rng.random_sample(size)*4.0 + 17.0
    i_ab = rng.random_sample(size)*4.0 + 17.0
    z_ab = rng.random_sample(size)*4.0 + 17.0
    y_ab = rng.random_sample(size)*4.0 + 17.0
    redshift = rng.random_sample(size)*2.0

    t0_mjd = mjd - rng.random_sample(size)*1000.0
    agn_tau = rng.random_sample(size)*1000.0 + 1000.0
    agnSeed = rng.random_integers(low=2, high=4000, size=size)
    agn_sfu = rng.random_sample(size)
    agn_sfg = rng.random_sample(size)
    agn_sfr = rng.random_sample(size)
    agn_sfi = rng.random_sample(size)
    agn_sfz = rng.random_sample(size)
    agn_sfy = rng.random_sample(size)

    rrStar = rng.random_sample(size)*np.radians(radius)
    thetaStar = rng.random_sample(size)*2.0*np.pi

    if deltaRA is None:
        raStar = centerRA + rrStar*np.cos(thetaStar)
    else:
        raStar = centerRA + np.radians(deltaRA)

    if deltaDec is None:
        decStar = centerDec + rrStar*np.sin(thetaStar)
    else:
        decStar = centerDec + np.radians(deltaDec)

    raStar = np.degrees(raStar)
    decStar = np.degrees(decStar)

    magnormStar = rng.random_sample(size)*4.0 + 17.0
    mudecl = rng.random_sample(size)*0.0001
    mura = rng.random_sample(size)*0.0001
    galacticAv = rng.random_sample(size)*0.05*3.1
    vrad = rng.random_sample(size)*1.0
    parallax = 0.00045+rng.random_sample(size)*0.00001
    period = rng.random_sample(size)*20.0
    amp = rng.random_sample(size)*5.0

    # write the data to the tables.
    for i in range(size):

        cmd = '''INSERT INTO galaxy_bulge VALUES (%i, %i, %f, %f, %f, %f, %f,
              '%s', %f, %f, %f, %i, '%s', %f, %f, %f, %f, %f, %f, %f, %f, %f, %f)''' % \
              (i, i, bra[i], bdec[i], ra[i], dec[i], magnorm_bulge[i], galaxy_seds[i%len(galaxy_seds)],
               a_b[i], b_b[i], pa_bulge[i], 4, 'CCM', av_b[i], rv_b[i], u_ab[i], g_ab[i],
               r_ab[i], i_ab[i], z_ab[i], y_ab[i], redshift[i], BulgeHalfLightRadius[i])

        c.execute(cmd)

        varParam = {'varMethodName': 'applyAgn',
                    'pars': {'agn_tau': round(agn_tau[i], 4), 't0_mjd': round(t0_mjd[i], 4),
                             'agn_sfu': round(agn_sfu[i], 4), 'agn_sfg': round(agn_sfg[i], 4),
                             'agn_sfr': round(agn_sfr[i], 4), 'agn_sfi': round(agn_sfi[i], 4),
                             'agn_sfz': round(agn_sfz[i], 4), 'agn_sfy': round(agn_sfy[i], 4),
                             'seed': int(agnSeed[i])}}

        paramStr = json.dumps(varParam)

        cmd = '''INSERT INTO galaxy VALUES (%i, %i, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f,
                                            '%s', '%s', '%s', '%s',
                                            %f, %f, %f, %i,
                                            %f, %f, %f, %i,
                                            '%s', %f, %f,
                                            '%s', %f, %f,
                                            %f, %f, %f, %f, %f, %f,
                                            %f, %f, %f)''' % \
              (i, i, ra[i], dec[i], bra[i], bdec[i], dra[i], ddec[i], agnra[i], agndec[i],
               magnorm_bulge[i], magnorm_disk[i], magnorm_agn[i],
               galaxy_seds[i%len(galaxy_seds)], galaxy_seds[i%len(galaxy_seds)], agn_sed,
               paramStr,
               a_b[i], b_b[i], pa_bulge[i], 4,
               a_d[i], b_d[i], pa_disk[i], 1,
               'CCM', av_b[i], rv_b[i],
               'CCM', av_d[i], rv_d[i],
               u_ab[i], g_ab[i], r_ab[i], i_ab[i], z_ab[i], y_ab[i], redshift[i],
               BulgeHalfLightRadius[i], DiskHalfLightRadius[i])
        c.execute(cmd)

        cmd = '''INSERT INTO galaxy_agn VALUES (%i, %i, %f, %f, %f, %f, %f, '%s', '%s',
              %f, %f, %f, %f, %f, %f, %f)''' % \
              (i, i, agnra[i], agndec[i], ra[i], dec[i],
               magnorm_agn[i], agn_sed, paramStr,
               u_ab[i], g_ab[i], r_ab[i], i_ab[i],
               z_ab[i], y_ab[i], redshift[i])

        c.execute(cmd)

        varParam = {'varMethodName': 'testVar', 'pars': {'period': period[i], 'amplitude': amp[i]}}
        paramStr = json.dumps(varParam)
        cmd = '''INSERT INTO StarAllForceseek VALUES (%i, %f, %f, %f, %f, %f, %f, %f, '%s', '%s', %f)''' %\
              (i, raStar[i], decStar[i], magnormStar[i], mudecl[i], mura[i],
               galacticAv[i], vrad[i], paramStr, star_seds[i%len(star_seds)], parallax[i])

        c.execute(cmd)

    conn.commit()
    conn.close()
    return obs_metadata
