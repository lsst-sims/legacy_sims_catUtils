"""
This is script will generate basic catalogs from the SSM tables on
fatboy in order to verify that the interface between SolarSystemObj
and fatboy is up to date.

This is not included as a unit test, because it depends on a connection
to fatboy and:

1) the connection to fatboy is still unreliable

2) we have not yet arranged for the Continuous Integration machines (e.g.
buildbot, Jenkins) to have a connection to fatboy.
"""

import numpy as np
import os
from lsst.utils import getPackageDir
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.baseCatalogModels import SolarSystemObj

class ssmBaseCatalog(InstanceCatalog):
    column_outputs = ['objid', 'raJ2000', 'decJ2000', 'sedFilename',
                      'velRa', 'velDec']

    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees,
                       'velRa': np.degrees, 'velDec': np.degrees}

if __name__ == "__main__":

    scratchDir = os.path.join(getPackageDir('sims_catUtils'),
                              'examples', 'scratch')

    catName = os.path.join(scratchDir, 'ssm_basic_catalog.txt')

    if os.path.exists(catName):
        os.unlink(catName)

    mjd = 50125.0

    obs = ObservationMetaData(unrefractedRA=25.0, unrefractedDec=-5.0,
                              mjd=mjd, boundType='circle',
                              boundLength=0.5)

    db = SolarSystemObj()

    cat = ssmBaseCatalog(db, obs_metadata=obs)

    cat.write_catalog(catName)

    typeList = [('id', np.int),
                ('ra', np.float),
                ('dec', np.float),
                ('name', str, 100),
                ('velRa', np.float),
                ('velDec', np.float)]

    dtype = np.dtype(typeList)

    data = np.genfromtxt(catName, dtype)
    for name in typeList:
        assert(len(data[name[0]])>0)

