"""
This script demonstrates how to use the stand-along findChipName method from astrometry
using an ObservationMetaData generated from the OpSim database
"""

import numpy
import os
from collections import OrderedDict
from lsst.sims.catalogs.generation.db import ObservationMetaData, \
                                             calcObsDefaults, getRotTelPos, \
                                             altAzToRaDec, Site
from lsst.sims.coordUtils import CameraCoords
from lsst.afw.cameraGeom.cameraConfig import CameraConfig
from lsst.afw.cameraGeom.cameraFactory import makeCameraFromPath
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels import OpSim3_61DBObject

def convertAmpNames(rawName):
   #This method takes the names of detectors as stored in a CameraConfig object
   #and converts them into the roots of the fits files storing amplifier information

   name = rawName.replace(':','')
   name = name.replace(',','')
   name = name.replace(' ','')
   name = name.replace('S','_S')
   return name

def makeLSSTcamera():
    filedir = os.path.join(os.getenv('OBS_LSSTSIM_DIR'), 'description/camera/')
    #this directory contains camera.py, which describes the lay out of
    #detectors on the LSST camera, as well as fits files that describe
    #the amplifiers on those detectors

    filename = os.path.join(filedir, 'camera.py')
    myCameraConfig = CameraConfig()
    myCameraConfig.load(filename)
    #converts the camera.py file into a CameraConfig object

    camera = makeCameraFromPath(myCameraConfig, filedir, convertAmpNames)
    #reads in the CameraConfig file as well as the fits files describing the
    #amplifiers and makes an afwCameraGeom.Camera object out of them

    return camera

epoch = 2000.0

#generate an ObservationMetaData object based on an actual OpSim pointing
obshistid = 88625744
radiusDegrees = 3.0
OpSimDB = CatalogDBObject.from_objid('opsim3_61')
obs_metadata = OpSimDB.getObservationMetaData(obshistid, radiusDegrees, makeCircBounds=True)

myCamCoords = CameraCoords()

camera = makeLSSTcamera()

#generate some random RA and Dec to find chips for
nsamples = 10
numpy.random.seed(32)
rr = numpy.radians(2.0)*numpy.random.sample(nsamples)
theta = 2.0*numpy.pi*numpy.random.sample(nsamples)
ra = numpy.radians(obs_metadata.unrefractedRA) + rr*numpy.cos(theta)
dec = numpy.radians(obs_metadata.unrefractedDec) + rr*numpy.sin(theta)

chipNames = myCamCoords.findChipName(ra=ra, dec=dec, epoch=epoch, camera=camera, obs_metadata=obs_metadata)

for (rr,dd,nn) in zip(ra,dec,chipNames):
    print rr,dd,nn
