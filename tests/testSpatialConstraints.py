from __future__ import print_function
import unittest
import numpy
import lsst
import lsst.utils.tests
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import haversine

# The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm

_testSpatialConstraints_is_connected = True

try:
    _test_db = bcm.GalaxyObj()
except:
    _testSpatialConstraints_is_connected = False


def setup_module(module):
    lsst.utils.tests.init()


class testCatalogBounds(unittest.TestCase):
    #@unittest.skipIf(not _testSpatialConstraints_is_connected,
    #                 "We are not connnected to fatboy")
    @unittest.skip('will fail until we have restored SSO tables on fatboy')
    def testCircleBounds(self):
        """Test Sql Server circular search region.
        exepectedFailure used despite expectation of success
        because the test depends on a network connection.
        """
        column_outputs = ['raJ2000', 'decJ2000']
        for objname, objcls in CatalogDBObject.registry.iteritems():
            if (not objcls.doRunTest or
                (objcls.testObservationMetaData is None) or
                (objcls.testObservationMetaData.bounds is None) or
                (objcls.testObservationMetaData.bounds.boundType != 'circle')):

                continue

            print("Running tests for", objname)
            obs_metadata = objcls.testObservationMetaData
            dbobj = objcls(verbose=False)
            result = dbobj.query_columns(column_outputs, obs_metadata=obs_metadata)

            # testObservationMetadata gives few enough results for one chunk
            try:
                result = next(result)
            except StopIteration:
                raise RuntimeError("No results for %s."%(objname))

            # confirm radius > distance from all points to center
            self.assertGreater(obs_metadata.bounds.radius + 1.e-4,
                               max(haversine(numpy.radians(obs_metadata.pointingRA),
                                             numpy.radians(obs_metadata.pointingDec),
                                             result['raJ2000'], result['decJ2000'])))

    #@unittest.skipIf(not _testSpatialConstraints_is_connected,
    #                 "We are not connected to fatboy")
    @unittest.skip('will fail until we have restored SSO tables on fatboy')
    def testBoxBounds(self):
        """Test Sql Server rectangular search region (ra/dec cuts).
        exepectedFailure used despite expectation of success
        because test depends on a network connection.
        """
        column_outputs = ['raJ2000', 'decJ2000']
        for objname, objcls in CatalogDBObject.registry.iteritems():
            if (not objcls.doRunTest or
                (objcls.testObservationMetaData is None) or
                (objcls.testObservationMetaData.bounds is None) or
                (objcls.testObservationMetaData.bounds.boundType != 'circle')):

                continue

            print("Running tests for", objname)
            circ_bounds = objcls.testObservationMetaData.bounds
            length = numpy.degrees(circ_bounds.radius)
            raCenter = numpy.degrees(circ_bounds.RA)+length
            decCenter = numpy.degrees(circ_bounds.DEC)+length
            obs_metadata = ObservationMetaData(boundType='box',
                                               pointingRA=raCenter,
                                               pointingDec=decCenter,
                                               boundLength=length,
                                               mjd=51000., bandpassName='i')
            dbobj = objcls(verbose=False)
            result = dbobj.query_columns(column_outputs, obs_metadata=obs_metadata)
            # testObservationMetadata gives few enough results for one chunk
            try:
                result = next(result)
            except StopIteration:
                raise RuntimeError("No results for %s."%(objname))

            self.assertLess(max(result['raJ2000']), numpy.radians(raCenter+length))
            self.assertGreater(min(result['raJ2000']), numpy.radians(raCenter-length))

            self.assertLess(max(result['decJ2000']), numpy.radians(decCenter+length))
            self.assertGreater(max(result['decJ2000']), numpy.radians(decCenter-length))


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
