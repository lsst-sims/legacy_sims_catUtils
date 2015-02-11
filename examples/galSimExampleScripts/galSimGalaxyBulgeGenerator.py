"""
This script shows how to use our GalSim interface to create FITS images of
galaxy bulges (or any sersic profile)
"""

import os
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catUtils.galSimInterface import *

#if you want to use the actual LSST camera
#from lsst.obs.lsstSim import LsstSimMapper

class testGalSimGalaxies(GalSimGalaxies):
    #only draw images for u and g bands (for speed)
    band_pass_names = ['u','g']

    #Note, we are not convolving with any PSF
    #see galSimStarGenerator.py


#select an OpSim pointing
obsMD = OpSim3_61DBObject()
obs_metadata = obsMD.getObservationMetaData(88625744, 0.05, makeCircBounds = True)

#grab a database of galaxies (in this case, galaxy bulges)
gals = CatalogDBObject.from_objid('galaxyBulge')

#now append a bunch of objects with 2D sersic profiles to our output file
galaxy_galSim = testGalSimGalaxies(gals, obs_metadata=obs_metadata)

#If you want to use the LSST camera, uncomment the line below.
#You can similarly assign any camera object you want here, as long
#as you do it before calling write_catalog()
#galaxy_galSim.camera = LsstSimMapper().camera


galaxy_galSim.write_catalog('galSim_bulge_example.txt', chunk_size=10000)
galaxy_galSim.write_images(nameRoot='bulge')
