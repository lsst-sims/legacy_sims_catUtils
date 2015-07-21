"""
This is a version of tutorial00.ipynb without the
running commentary
"""

from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels import *

myGalaxyDB = CatalogDBObject.from_objid('galaxyTiled')
myStarDB = CatalogDBObject.from_objid('allstars')

from lsst.sims.utils import ObservationMetaData

obs_metadata = ObservationMetaData(unrefractedRA = 220.0,
                                   unrefractedDec = 19.0,
                                   boundType = 'circle',
                                   boundLength = 0.2,
                                   mjd = 52000.0)


from lsst.sims.catUtils.exampleCatalogDefinitions import RefCatalogGalaxyBase, \
                                                  RefCatalogStarBase

myStarCat = RefCatalogStarBase(myStarDB, obs_metadata=obs_metadata)
myStarCat.write_catalog('star_example.txt')

myGalaxyCat = RefCatalogGalaxyBase(myGalaxyDB, obs_metadata=obs_metadata)
myGalaxyCat.write_catalog('galaxy_example.txt')


squareObsMetadata = ObservationMetaData(unrefractedRA = 220.0,
                                       unrefractedDec = 19.0,
                                       boundType = 'box',
                                       boundLength = 0.3,
                                       mjd = 52000.0)

myStarCat = RefCatalogStarBase(myStarDB, obs_metadata=squareObsMetadata)
myStarCat.write_catalog('star_example_square.txt')
