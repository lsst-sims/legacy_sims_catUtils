from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.utils import ObservationMetaData
#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catUtils.exampleCatalogDefinitions import RefCatalogGalaxyBase

if __name__ == '__main__':

    obs_metadata = ObservationMetaData(boundType='circle', pointingRA=0.0, pointingDec=0.0,
                                       boundLength=0.01, mjd=57388.0)
    dbobj = CatalogDBObject.from_objid('galaxyBase')
    t = dbobj.getCatalog('galaxy_photometry_cat', obs_metadata=obs_metadata)
    filename = 'galaxy_photometry_test.dat'
    t.write_catalog(filename, chunk_size=10000)

