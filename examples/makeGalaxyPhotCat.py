from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catUtils.exampleCatalogDefinitions import RefCatalogGalaxyBase

if __name__ == '__main__':

    obs_metadata = ObservationMetaData(boundType='circle', unrefractedRA=0.0, unrefractedDec=0.0,
                                       boundLength=0.01)
    dbobj = DBObject.from_objid('galaxyBase')
    t = dbobj.getCatalog('galaxy_photometry_cat', obs_metadata=obs_metadata)
    filename = 'galaxy_photometry_test.dat'
    t.write_catalog(filename, chunk_size=10000)

