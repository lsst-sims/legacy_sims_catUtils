from __future__ import with_statement
from __future__ import print_function
import os
import sys
import traceback
import unittest
import lsst.utils.tests

from lsst.sims.catUtils.utils import failedOnFatboy

import numpy as np
from lsst.sims.utils.CodeUtilities import sims_clean_up
# Observation metadata modules
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
# To access opsim sqlite database
from lsst.utils import getPackageDir
# photometric parameters (exptime lives here for dmag calculation)
from lsst.sims.photUtils import PhotometricParameters
# SSM catalog modules
from lsst.sims.catUtils.baseCatalogModels import SolarSystemObj
from lsst.sims.catUtils.mixins import PhotometrySSM, AstrometrySSM, CameraCoords, ObsMetadataBase
from lsst.sims.catalogs.definitions import InstanceCatalog
# For camera.
from lsst.obs.lsstSim import LsstSimMapper

import time


def setup_module(module):
    lsst.utils.tests.init()


def dtime(time_prev):
    return (time.time() - time_prev, time.time())


def reassure():
    print('\ntestObsCat failed to connect to fatboy')
    print('Sometimes that happens.  Do not worry.')

# Build sso instance class
basic_columns = ['objid', 'expMJD', 'raJ2000', 'decJ2000', 'velRa', 'velDec',
                 'skyVelocity', 'dist', 'dmagTrailing', 'dmagDetection',
                 'sedFilename', 'magFilter', 'magSNR', 'visibility',
                 'seeing', 'bandpass', 'visitExpTime', 'm5']


class ssmCat(InstanceCatalog, PhotometrySSM, AstrometrySSM, ObsMetadataBase, CameraCoords):
    catalog_type = __file__ + 'ssm_cat'

    column_outputs = basic_columns
    cannot_be_null = ['visibility']
    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees,
                       'velRa': np.degrees, 'velDec': np.degrees}
    default_formats = {'f': '%.13f'}


class ssmCatCamera(ssmCat):
    catalog_type = __file__ + 'ssm_cat_camera'

    column_outputs = basic_columns + ['chipName']
    camera = LsstSimMapper().camera
    cannot_be_null = ['visibility', 'chipName']
    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees,
                       'velRa': np.degrees, 'velDec': np.degrees}
    default_formats = {'f': '%.13f'}

######


class createSSMSourceCatalogsTest(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()

    def test_ssm_catalog_creation(self):

        t = time.time()
        # Fake opsim data.
        database = os.path.join(getPackageDir('SIMS_DATA'), 'OpSimData/opsimblitz1_1133_sqlite.db')
        generator = ObservationMetaDataGenerator(database=database, driver='sqlite')

        night = 20
        query = 'select min(expMJD), max(expMJD) from summary where night=%d' % (night)
        res = generator.opsimdb.execute_arbitrary(query)
        expMJD_min = res[0][0]
        expMJD_max = res[0][1]

        obsMetaDataResults = generator.getObservationMetaData(expMJD=(expMJD_min, expMJD_max),
                                                              limit=3, boundLength=2.2)

        dt, t = dtime(t)
        print('To query opsim database: %f seconds' % (dt))

        write_header = True
        write_mode = 'w'

        try:
            ssmObj = SolarSystemObj()

            for obsMeta in obsMetaDataResults:
                # But moving objects databases are not currently complete for all years.
                # Push forward to night=747.
                # (note that we need the phosim dictionary as well)
                newMJD = obsMeta.mjd.TAI + (747 - 20)
                phoSimMetaDict = {'exptime': [30]}
                obs = ObservationMetaData(mjd=newMJD,
                                          pointingRA=obsMeta.pointingRA,
                                          pointingDec=obsMeta.pointingDec,
                                          bandpassName=obsMeta.bandpass,
                                          rotSkyPos=obsMeta.rotSkyPos,
                                          m5=obsMeta.m5[obsMeta.bandpass],
                                          seeing=obsMeta.seeing[obsMeta.bandpass],
                                          boundLength=obsMeta.boundLength,
                                          boundType=obsMeta.boundType)

                obs.phoSimMetaData = phoSimMetaDict

                mySsmDb = ssmCatCamera(ssmObj, obs_metadata = obs)
                photParams = PhotometricParameters(exptime = obs.phoSimMetaData['exptime'][0],
                                                   nexp=1, bandpass=obs.bandpass)
                mySsmDb.photParams = photParams

                try:
                    with lsst.utils.tests.getTempFilePath('.txt') as output_cat:
                        mySsmDb.write_catalog(output_cat, write_header=write_header, write_mode=write_mode)

                        # verify that we did not write an empty catalog
                        with open(output_cat, 'r') as input_file:
                            lines = input_file.readlines()
                        self.assertGreater(len(lines), 1)
                except:
                    # This is because the solar system object 'tables'
                    # don't actually connect to tables on fatboy; they just
                    # call methods stored on fatboy.  Therefore, the connection
                    # failure will not be noticed until this part of the test
                    msg = sys.exc_info()[1].args[0]
                    if 'DB-Lib error' in msg:
                        reassure()
                        continue
                    else:
                        raise

                write_mode = 'a'
                write_header = False

                dt, t = dtime(t)
                print('To query solar system objects: %f seconds (obs MJD time %f)' % (dt, obs.mjd.TAI))

        except:
            trace = traceback.extract_tb(sys.exc_info()[2], limit=20)
            msg = sys.exc_info()[1].args[0]
            if 'Failed to connect' in msg or failedOnFatboy(trace):
                # if the exception was because of a failed connection
                # to fatboy, ignore it.
                reassure()

                pass
            else:
                raise


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
