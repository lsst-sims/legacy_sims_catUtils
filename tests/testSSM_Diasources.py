from __future__ import with_statement
import os, sys
import traceback
import unittest
from copy import deepcopy
import lsst.utils.tests as utilsTests

import numpy as np
# Observation metadata modules
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
# To access opsim sqlite database
from lsst.utils import getPackageDir
# photometric parameters (exptime lives here for dmag calculation)
from lsst.sims.photUtils import PhotometricParameters
# SSM catalog modules
from lsst.sims.catUtils.baseCatalogModels import SolarSystemObj, CometObj, MBAObj, NEOObj, MiscSolarSystemObj
from lsst.sims.catUtils.mixins import PhotometrySSM, AstrometrySSM, CameraCoords, ObsMetadataBase
from lsst.sims.catalogs.measures.instance import InstanceCatalog
# For camera.
from lsst.obs.lsstSim import LsstSimMapper

import time
def dtime(time_prev):
    return (time.time() - time_prev, time.time())

def failedOnFatboy(tracebackList):
    """
    Accepts a list generated by traceback.extract_tb; determines if the last
    point in the sims code in the traceback is _connect_to_engine (from
    sims_catalogs_generation/../db/dbConnection.py), in which case, the failure
    was probably due to fatboy connectivity.
    """
    if not isinstance(tracebackList, list):
        return False

    lastSimsDex = -1
    for ix, item in enumerate(tracebackList):
        if not isinstance(item, tuple):
            return False

        if 'sims' in item[0]:
            lastSimsDex = ix

    if lastSimsDex<0:
        return False

    if '_connect_to_engine' in tracebackList[lastSimsDex][2]:
        return True

    return False


def reassure():
    print '\ntestObsCat failed to connect to fatboy'
    print 'Sometimes that happens.  Do not worry.'

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

class createSSMSourceCatalogsTest(unittest.TestCase):

    def test_ssm_catalog_creation(self):

        t = time.time()
        # Fake opsim data.
        database = os.path.join(getPackageDir('SIMS_DATA'), 'OpSimData/opsimblitz1_1133_sqlite.db')
        generator = ObservationMetaDataGenerator(database=database, driver='sqlite')

        night = 20
        query = 'select min(expMJD), max(expMJD) from summary where night=%d' %(night)
        res = generator.opsimdb.execute_arbitrary(query)
        expMJD_min = res[0][0]
        expMJD_max = res[0][1]

        obsMetaDataResults = generator.getObservationMetaData(expMJD=(expMJD_min, expMJD_max), limit=3, boundLength=2.2)

        dt, t = dtime(t)
        print 'To query opsim database: %f seconds' %(dt)

        write_header = True
        write_mode = 'w'

        try:
            #ssmObj = NEOObj()
            ssmObj = SolarSystemObj()

            output_cat = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace', 'catsim_ssm_test')
            if os.path.exists(output_cat):
                os.unlink(output_cat)

            for obsMeta in obsMetaDataResults:
                # But moving objects databases are not currently complete for all years. Push forward to night=747.
                # (note that we need the phosim dictionary as well)
                newMJD = obsMeta.mjd.TAI + (747 - 20)
                phoSimMetaDict = {'exptime': [30]}
                obs = ObservationMetaData(phoSimMetaData = phoSimMetaDict, mjd=newMJD,
                                          pointingRA=obsMeta.pointingRA, pointingDec=obsMeta.pointingDec,
                                          bandpassName=obsMeta.bandpass, rotSkyPos=obsMeta.rotSkyPos,
                                          m5=obsMeta.m5[obsMeta.bandpass], seeing=obsMeta.seeing[obsMeta.bandpass],
                                          boundLength=obsMeta.boundLength, boundType=obsMeta.boundType)
                mySsmDb = ssmCatCamera(ssmObj, obs_metadata = obs)
                #mySsmDb = ssmCat(ssmObj, obs_metadata = obs)
                photParams = PhotometricParameters(exptime = obs.phoSimMetaData['exptime'][0], nexp=1, bandpass=obs.bandpass)
                mySsmDb.photParams = photParams

                try:
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
                print 'To query solar system objects: %f seconds (obs MJD time %f)' %(dt, obs.mjd.TAI)

                if os.path.exists(output_cat):
                    os.unlink(output_cat)

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


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(createSSMSourceCatalogsTest)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
