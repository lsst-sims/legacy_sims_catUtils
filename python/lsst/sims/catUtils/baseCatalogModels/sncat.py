#!/usr/bin/env python

from cStringIO import StringIO
import sys
import os
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.generation.db import ObservationMetaData
from lsst.sims.photUtils import Sed
import lsst.sims.photUtils.Bandpass as Bandpass
import sncosmo
from astropy.units import Unit
import astropy.cosmology as cosmology 
#from astropy.cosmology import Planck13 as cosmo
from snObject import SNObject
from lsst.sims.photUtils.CosmologyObject import CosmologyWrapper 
import sqliteutils as sq
import sqlite3
wavelenstep = 0.1



# class SNIaCatalog (object):
class SNIaCatalog (InstanceCatalog):
    """
    Supernova Type Ia in the catalog are characterized by the  following
    attributes
    position (ra, dec, redshift),
    velocity wrt host galaxy,
    the supernova model (eg. SALT2)
    and parameters of the supernova model that predict the SED.
    """
    column_outputs = ['snid', 'snra', 'sndec', 'z', 't0', 'c', 'x1',
                      'x0', 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z',
                      'mag_y']
    override_formats = {'snra': '%8e', 'sndec': '%8e', 'c': '%8e',
                        'x0': '%8e'}
    cannot_be_null = ['x0']
# column_outputs=['raJ2000','decJ2000','snid','z','snra', 'sndec',\
# 'mass_stellar', 'c', 'x1', 't0', "x0"]
    surveyoffset = 570000.0
    SN_thresh = 100.0
    def usedlsstbands(self, loadsncosmo=True, loadcatsim=True):

        bandPassList = self.obs_metadata.bandpass
        banddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')
        lsstbands = []
        lsstbp = {}

        # wavelenstep = 0.1
        for band in bandPassList:
            # setup sncosmo bandpasses
            bandfname = banddir + "/total_" + band + '.dat'
            if loadsncosmo:
                # Usually the next two lines can be merged,
                # but there is an astropy bug currently.
                numpyband = np.loadtxt(bandfname)
                sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
                                               trans=numpyband[:, 1],
                                               wave_unit=Unit('nm'),
                                               name='LSST'+band)
                sncosmo.registry.register(sncosmoband, force=True)
            if loadcatsim:
                # Now load LSST bandpasses for catsim
                lsstbp[band] = Bandpass()
                lsstbp[band].readThroughput(bandfname,
                                            wavelen_step=wavelenstep)
                lsstbands.append(lsstbp[band])
        plot = False
        if plot:
            filterfigs, filterax = plt.subplots()
            for band in bandPassList:
                b = sncosmo.get_bandpass('LSST' + band)
                filterax.plot(b.wave, b.trans, '-k', lw=2.0)
                filterax.set_xlabel(r'$\lambda$, ($\AA$)')
                filterax.set_ylabel(r'transmission')
            plt.show()

        if loadcatsim:
            return lsstbands
        else:
            return None

    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        # rethink
        return self.column_by_name('id')

    @property
    def numobjs(self):
        return len(self.column_by_name('id'))

    @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
    def get_angularCoordinates(self):
        _snra, _sndec, _z = self.column_by_name('raJ2000'), \
            self.column_by_name('decJ2000'), \
            self.column_by_name('redshift')
        _sndec += np.zeros(self.numobjs)
        _snra += np.zeros(self.numobjs)
        _vra = np.zeros(self.numobjs)
        _vdec = np.zeros(self.numobjs)
        _vr = np.zeros(self.numobjs)
        return ([_snra, _sndec, _z, _vra, _vdec, _vr])

    @compound('c', 'x1', 'x0', 't0', 'mag_u',
              'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y')
    def get_snparams(self):
        lsstbands = self.usedlsstbands()
        ra, dec = self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')
        SNmodel = SNObject()
        hundredyear = 100*365.0
        vals = np.zeros(shape=(self.numobjs, 10))
        _z, _id = self.column_by_name('redshift'), self.column_by_name('snid')
        #distmods = cosmo.get_cosmologicalDistanceModulus()
        bad = np.nan
        for i, v in enumerate(vals):
            np.random.seed(_id[i])
            v[-1] = np.random.uniform(-hundredyear / 2.0 +
                                      self.surveyoffset,
                                      hundredyear / 2.0 +
                                      self.surveyoffset)
            if np.abs(v[-1] - self.obs_metadata.mjd) > self.SN_thresh:
                # v = np.array([np.nan, np.nan, np.nan, np.nan])
                v[-1] = bad
                v[0] = bad
                v[1] = bad
                v[2] = bad
                continue
            v[0] = np.random.normal(0., 0.3)
            v[1] = np.random.normal(0., 3.0)
            mabs = np.random.normal(-19.3, 0.3)
            SNmodel.ra = ra[i]
            SNmodel.dec = dec[i]
            SNmodel.mwebvfrommaps()
            SNmodel.set(z=_z[i], c=v[0], x1=v[1], t0=v[-1])
            # rather than use the SNCosmo function below which uses astropy to calculate
            # distanceModulus, we will use photUtils CosmologyWrapper for consistency
            # SNmodel.set_source_peakabsmag(mabs, 'bessellb', 'ab', cosmo=cosmo)
            mag = mabs + cosmo.cosmology.distanceModulus(_z[i])
            #cosmo.get_cosmologicalDistanceModulus(_z[i])
            SNmodel.source.set_peakmag(mag, band='bessellb', magsys='ab')
            v[2] = SNmodel.get('x0')
            v[3] = v[-1]
            v[4:] = SNmodel.lsstbandmags(lsstbands=lsstbands,
                                         time=self.obs_metadata.mjd)
            # print self.obs_metadata.mjd
        # print self.obs_metadata.bandpass

        return ([vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3],
                vals[:, 4], vals[:, 5], vals[:, 6], vals[:, 7],
                vals[:, 8], vals[:, 9]])

if __name__ == "__main__":

    import lsst.sims.catUtils.baseCatalogModels as bcm
    # import timeit
    print bcm.__file__
    from lsst.sims.catalogs.generation.db import ObservationMetaData


    # Define the cosmology for the catalog
    cosmo = CosmologyWrapper() 
    # You can set the cosmology w0waCDM FLRW cosmologies
    cosmo.setCosmology(Om0=0.25, Ok0=None, H0=73.0)

    galDB = CatalogDBObject.from_objid('galaxyTiled')
    def file2lst(fname, i, mjd):
        d = np.loadtxt(fname, delimiter=',')
        l = list()
        for i, row in enumerate(d):
            obsid = 'obshist' + str(i)
            lst = [obsid] + [mjd] + row.tolist()
            l.append(lst)
        return l
    # main example :  create all the 'observations' separated by a day, 

    myMJDS = [570123.15 + 3.*i for i in range(20)]

    createcat = True
    #prepare a new sncat table:
    #delete catalog if it exists manually
    connection = sqlite3.connect('../out/sncat.db')
    curs = connection.cursor()
    curs.execute('CREATE TABLE if not exists mysncat (id TEXT, mjd FLOAT, snid INT, snra FLOAT, sndec FLOAT, z FLOAT, t0 FLOAT, c FLOAT, x1 FLOAT, x0 FLOAT, mag_u FLOAT, mag_g FLOAT, mag_r FLOAT, mag_i FLOAT, mag_z FLOAT, mag_y FLOAT)')
    for i, myMJD in enumerate(myMJDS):
        if createcat:
            myObsMD = ObservationMetaData(boundType='circle',
                                          unrefractedRA=5.0,
                                          unrefractedDec=15.0,
                                          boundLength=0.15,
                                          bandpassName=['u', 'g', 'r', 'i',
                                                        'z', 'y'],
                                          mjd=myMJD)
            catalog = SNIaCatalog(db_obj=galDB,
                                  obs_metadata=myObsMD)
            print "====================================="
            print i, type(catalog.usedlsstbands()) , catalog.obs_metadata.mjd
            print "====================================="
            fname = "../out/SNIaCat_" + str(i) + ".txt"
            catalog.write_catalog(fname)
        fname = "../out/SNIaCat_" + str(i) + ".txt"
        l = file2lst(fname, i, mjd=myMJD)
        recs = sq.array2dbrecords(l)
        exec_str = sq.insertfromdata(tablename='mysncat', records=recs, multiple=True)
        curs.executemany(exec_str, recs)
    connection.commit()
    curs.execute('SELECT * FROM mysncat')  
    print 'In Database: ', curs.fetchall()
    connection.close()


