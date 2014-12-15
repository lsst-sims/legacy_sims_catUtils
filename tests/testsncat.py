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
from sncat import SNIaCatalog
from lsst.sims.photUtils.CosmologyObject import CosmologyWrapper 
#import sqliteutils as sq
import utils_for_test as sq
import sqlite3
wavelenstep = 0.1

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
connection = sqlite3.connect('testData/sncat.db')
# connection = sqlite3.connect('testData/STDsncat.db')
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
        fname = "testData/SNIaCat_" + str(i) + ".txt"
        catalog.write_catalog(fname)
    fname = "testData/SNIaCat_" + str(i) + ".txt"
    l = file2lst(fname, i, mjd=myMJD)
    recs = sq.array2dbrecords(l)
    exec_str = sq.insertfromdata(tablename='mysncat', records=recs, multiple=True)
    curs.executemany(exec_str, recs)
connection.commit()
curs.execute('SELECT * FROM mysncat')  
print 'In Database: ', curs.fetchall()
connection.close()


