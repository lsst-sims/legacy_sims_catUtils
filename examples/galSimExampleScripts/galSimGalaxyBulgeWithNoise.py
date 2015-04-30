"""
This script shows how incorporate noise in images of galaxies
"""

import os
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import GalaxyBulgeObj, OpSim3_61DBObject
from lsst.sims.catUtils.galSimInterface import GalSimGalaxies, ExampleCCDNoise, \
                                               SNRdocumentPSF

#if you want to use the actual LSST camera
#from lsst.obs.lsstSim import LsstSimMapper

class testGalSimGalaxiesNoiseless(GalSimGalaxies):
    #only draw images for u and g bands (for speed)
    bandpass_names = ['u','g']

    #If you want to use the LSST camera, uncomment the line below.
    #You can similarly assign any camera object you want here
    #camera = LsstSimMapper().camera

    PSF = SNRdocumentPSF()

class testGalSimGalaxiesNoisy(testGalSimGalaxiesNoiseless):

    #defined in galSimInterface/galSimUtilities.py
    noise_and_background = ExampleCCDNoise(99)

#select an OpSim pointing
obsMD = OpSim3_61DBObject()
obs_metadata = obsMD.getObservationMetaData(88625744, 0.05, makeCircBounds = True)

#grab a database of galaxies (in this case, galaxy bulges)
gals = CatalogDBObject.from_objid('galaxyBulge')

#now append a bunch of objects with 2D sersic profiles to our output file
gal_noiseless = testGalSimGalaxiesNoiseless(gals, obs_metadata=obs_metadata)

gal_noiseless.write_catalog('galSim_NoiselessGalaxies_example.txt', chunk_size=10000)
gal_noiseless.write_images(nameRoot='noiselessGalaxies')

gal_noisy = testGalSimGalaxiesNoisy(gals, obs_metadata=obs_metadata)
gal_noisy.write_catalog('galSim_NoisyGalaxies_example.txt', chunk_size=10000)
gal_noisy.write_images(nameRoot='noisyGalaxies')
