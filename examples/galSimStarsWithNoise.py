"""
This script illustrates how to add noise to a series of FITS images using stars
"""

import os
import galsim
from lsst.sims.catalogs.generation.db import radiansToArcsec
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catUtils.galSimInterface import *

#if you want to use the actual LSST camera
#from lsst.obs.lsstSim import LsstSimMapper

class testGalSimStarsNoiseless(GalSimStars):
    #only draw images for u and g bands (for speed)
    band_pass_names = ['u','g']

    PSF = ExampleOpticalPSF()

class testGalSimStarsWithNoise(testGalSimStarsNoiseless):

    noise = ExampleCCDNoise(seed=99)

#select an OpSim pointing
obsMD = OpSim3_61DBObject()
obs_metadata = obsMD.getObservationMetaData(88625744, 0.1, makeCircBounds = True)

#grab a database of galaxies (in this case, galaxy bulges)
stars = CatalogDBObject.from_objid('allstars')

#now append a bunch of objects with 2D sersic profiles to our output file
stars_noiseless = testGalSimStarsNoiseless(stars, obs_metadata=obs_metadata)

#If you want to use the LSST camera, uncomment the line below.
#You can similarly assign any camera object you want here, as long
#as you do it before calling write_catalog()
#stars_noiseless.camera = LsstSimMapper().camera

stars_noiseless.write_catalog('galSim_NoiselessStars_example.txt', chunk_size=10000)
stars_noiseless.write_images(nameRoot='noiselessStars')

stars_noisy = testGalSimStarsWithNoise(stars, obs_metadata=obs_metadata)
stars_noisy.write_catalog('galSim_NoisyStars_example.txt', chunk_size=10000)
stars_noisy.add_noise()
stars_noisy.write_images(nameRoot='noisyStars')
