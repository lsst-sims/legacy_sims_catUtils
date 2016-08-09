"""
This is script will generate basic catalogs from the SSM tables on
fatboy in order to verify that the interface between SolarSystemObj
and all of its daughter classes and fatboy is up to date.

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
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.baseCatalogModels import SolarSystemObj, CometObj, \
                                                 MBAObj, NEOObj, \
                                                 MiscSolarSystemObj

class ssmBaseCatalog(InstanceCatalog):
    column_outputs = ['objid', 'raJ2000', 'decJ2000', 'sedFilename',
                      'velRa', 'velDec']

    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees,
                       'velRa': np.degrees, 'velDec': np.degrees}

if __name__ == "__main__":

    scratchDir = os.path.join(getPackageDir('sims_catUtils'),
                              'examples', 'scratch')

    catDict = {}
    catDict['ssm_basic_catalog.txt'] = SolarSystemObj
    catDict['ssm_comet_catalog.txt'] = CometObj
    catDict['ssm_neo_catalog.txt'] = NEOObj
    catDict['ssm_mba_catalog.txt'] = MBAObj
    catDict['ssm_misc_catalog.txt'] = MiscSolarSystemObj

    mjd = 50125.0

    for name in catDict:

        catName = os.path.join(scratchDir, name)

        if os.path.exists(catName):
            os.unlink(catName)


        obs = ObservationMetaData(pointingRA=25.0, pointingDec=-5.0,
                                  mjd=mjd, boundType='circle',
                                  boundLength=0.5)

        db = catDict[name]()

        cat = ssmBaseCatalog(db, obs_metadata=obs)

        cat.write_catalog(catName)

        with open(catName,'r') as readFile:
            lines = readFile.readlines()
            if len(lines) <= 1:
                raise RuntimeError('%s is empty' % catName)


