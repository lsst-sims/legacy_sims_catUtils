from __future__ import print_function

import numpy as np
# Observation metadata modules
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.photUtils import PhotometricParameters
# SSM catalog modules
from lsst.sims.catUtils.baseCatalogModels import SolarSystemObj, CometObj, MBAObj, NEOObj, MiscSolarSystemObj
from lsst.sims.catUtils.mixins import PhotometrySSM, AstrometrySSM, CameraCoords, ObsMetadataBase
from lsst.sims.catalogs.definitions import InstanceCatalog
# For camera.
from lsst.obs.lsstSim import LsstSimMapper

import time
def dtime(time_prev):
    return (time.time() - time_prev, time.time())

# Build sso instance class

basic_columns = ['objid', 'expMJD', 'raJ2000', 'decJ2000', 'velRa', 'velDec', 'skyVelocity', 'dist', 'dmagTrailing', 'dmagDetection',
                 'sedFilename', 'magFilter', 'magSNR', 'visibility', 'seeing', 'bandpass', 'visitExpTime', 'm5']

class ssmCat(InstanceCatalog, PhotometrySSM, AstrometrySSM, ObsMetadataBase, CameraCoords):
    column_outputs = basic_columns
    cannot_be_null = ['visibility']
    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees,
                       'velRa': np.degrees, 'velDec': np.degrees}
    default_formats = {'f':'%.13f'}


class ssmCatCamera(ssmCat):
    column_outputs = basic_columns + ['chipName']
    camera = LsstSimMapper().camera
    cannot_be_null = ['visibility', 'chipName']
    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees,
                       'velRa': np.degrees, 'velDec': np.degrees}
    default_formats = {'f':'%.13f'}

######

t = time.time()
# Get opsim data.
opsdb = '/Users/lynnej/opsim/db/minion_1016_newsky.db'
generator = ObservationMetaDataGenerator(database=opsdb, driver='sqlite')

night = 203
query = 'select min(expMJD), max(expMJD) from summary where night=%d and filter="r"' %(night)
res = generator.opsimdb.execute_arbitrary(query)
expMJD_min = res[0][0]
expMJD_max = res[0][1]

# Test image (deep, r band, near ecliptic)
obsMetaDataResults = generator.getObservationMetaData(expMJD=(expMJD_min, expMJD_max),  boundLength=2.2)

dt, t = dtime(t)
print('To query opsim database: %f seconds' %(dt))

write_header = True
write_mode = 'w'

#ssmObj = NEOObj()
ssmObj = SolarSystemObj()

for obs in obsMetaDataResults:
    #print obs.mjd, obs.unrefractedRA, obs.unrefractedDec, obs.bandpass, obs.boundType, obs.boundLength

    mySsmDb = ssmCatCamera(ssmObj, obs_metadata = obs)
    photParams = PhotometricParameters(exptime = obs.OpsimMetaData['visitExpTime'], nexp=1, bandpass=obs.bandpass)
    mySsmDb.photParams = photParams

    mySsmDb.write_catalog('catsim_n%d' %(night), write_header=write_header, write_mode=write_mode)
    write_mode = 'a'
    write_header = False

    dt, t = dtime(t)
    print('To query solar system objects: %f seconds (obs time %f)' %(dt, obs.mjd.TAI))
