import os
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catUtils.galSimInterface import *

class testGalSimGalaxies(GalSimGalaxies):
    band_pass_names = ['u','g']

#if you want to use the actual LSST camera
#from lsst.obs.lsstSim import LsstSimMapper

#select an OpSim pointing
obsMD = OpSim3_61DBObject()
obs_metadata = obsMD.getObservationMetaData(88625744, 0.01, makeCircBounds = True)

#grab a database of galaxies (in this case, galaxy bulges)
gals = CatalogDBObject.from_objid('galaxyBulge')

#now append a bunch of objects with 2D sersic profiles to our output file
galaxy_galSim = testGalSimGalaxies(gals, obs_metadata=obs_metadata)

#If you want to use the LSST camera, uncomment the line below.
#You can similarly assign any camera object you want here, as long
#as you do it immediately after instantiating GalSimGalaxies
#galaxy_galSim.camera = LsstSimMapper().camera


galaxy_galSim.write_catalog('galSim_bulge_example.txt', chunk_size=100)
galaxy_galSim.write_images(nameRoot='bulge')
