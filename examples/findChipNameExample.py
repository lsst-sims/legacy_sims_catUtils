"""
This script demonstrates how to use the stand-along findChipName method from astrometry
using an ObservationMetaData generated from the OpSim database
"""

import numpy
import os
from collections import OrderedDict
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.utils import pupilCoordsFromRaDec, observedFromICRS
from lsst.sims.coordUtils import chipNameFromRaDec
from lsst.sims.catUtils.baseCatalogModels import OpSim3_61DBObject
from lsst.obs.lsstSim import LsstSimMapper

mapper = LsstSimMapper()
camera = mapper.camera

epoch = 2000.0

#generate an ObservationMetaData object based on an actual OpSim pointing
obshistid = 88625744
radiusDegrees = 3.0
OpSimDB = OpSim3_61DBObject()
obs_metadata = OpSimDB.getObservationMetaData(obshistid, radiusDegrees, makeCircBounds=True)

#generate some random RA and Dec to find chips for
nsamples = 10
numpy.random.seed(32)
rr = numpy.radians(2.0)*numpy.random.sample(nsamples)
theta = 2.0*numpy.pi*numpy.random.sample(nsamples)
ra = numpy.degrees(numpy.radians(obs_metadata.pointingRA) + rr*numpy.cos(theta))
dec = numpy.degrees(numpy.radians(obs_metadata.pointingDec) + rr*numpy.sin(theta))

#need to correct coordinates for precession, nutation, and aberration
ra, dec = observedFromICRS(ra, dec, obs_metadata=obs_metadata, epoch=epoch)

xx, yy = pupilCoordsFromRaDec(ra, dec, obs_metadata=obs_metadata, epoch=epoch)

chipNames = chipNameFromRaDec(ra, dec, epoch=epoch, camera=camera, obs_metadata=obs_metadata)

for (rr,dd,x,y,nn) in zip(ra,dec,xx,yy,chipNames):
    print rr,dd,arcsecFromRadians(x),arcsecFromRadians(y),nn
