import unittest
import numpy
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData

#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm

def distanceToPoint(lon1, lat1, lon2, lat2):
    nx1 = numpy.cos(lat1) * numpy.cos(lon1)
    ny1 = numpy.cos(lat1) * numpy.sin(lon1)
    nz1 = numpy.sin(lat1)
    nx2 = numpy.cos(lat2) * numpy.cos(lon2)
    ny2 = numpy.cos(lat2) * numpy.sin(lon2)
    nz2 = numpy.sin(lat2)
    return 2 * numpy.arcsin(
                             numpy.sqrt((nx1 - nx2) * (nx1 - nx2)
                                      + (ny1 - ny2) * (ny1 - ny2)
                                      + (nz1 - nz2) * (nz1 - nz2)
                                        ) / 2.
                                        )
                        

vDistanceBetweenPoints = numpy.vectorize(distanceToPoint)

class testCatalogBounds(unittest.TestCase):
    @unittest.expectedFailure
    def testCircleBounds(self):
        """Test Sql Server circular search region.
        exepectedFailure used despite expectation of success
        because the test depends on a network connection.
        """
        column_outputs = ['raJ2000', 'decJ2000']
        for objname, objcls in CatalogDBObject.registry.iteritems():
            if not objcls.doRunTest \
            or (objcls.testObservationMetaData is None) \
            or (objcls.testObservationMetaData.bounds is None) \
            or (objcls.testObservationMetaData.bounds.boundType != 'circle'):
                continue
            print "Running tests for", objname
            obs_metadata = objcls.testObservationMetaData
            for i in range(5):
                try:
                    dbobj = objcls(verbose=False)
                    result = dbobj.query_columns(column_outputs, obs_metadata=obs_metadata)
                    break
                except:
                    print "network failure. Retrying."
            #testObservationMetadata gives few enough results for one chunk
            try:
                result = result.next()
            except StopIteration:
                raise RuntimeError("No results for %s."%(objname))

            #confirm radius > distance from all points to center
            self.assertGreater(obs_metadata.bounds.radius + 1.e-4,
                           max(vDistanceBetweenPoints(obs_metadata.unrefractedRA,
                                                      obs_metadata.unrefractedDec,
                                                      result['raJ2000'], result['decJ2000'])))

    @unittest.expectedFailure
    def testBoxBounds(self):
        """Test Sql Server rectangular search region (ra/dec cuts).
        exepectedFailure used despite expectation of success
        because test depends on a network connection.
        """
        column_outputs = ['raJ2000', 'decJ2000']
        for objname, objcls in CatalogDBObject.registry.iteritems():
            if not objcls.doRunTest \
            or (objcls.testObservationMetaData is None) \
            or (objcls.testObservationMetaData.bounds is None) \
            or (objcls.testObservationMetaData.bounds.boundType != 'circle'):
                continue
            print "Running tests for", objname
            circ_bounds = objcls.testObservationMetaData.bounds
            length = numpy.degrees(circ_bounds.radius)
            raCenter = numpy.degrees(circ_bounds.RA)+length
            decCenter = numpy.degrees(circ_bounds.DEC)+length
            obs_metadata = ObservationMetaData(boundType='box',unrefractedRA=raCenter,unrefractedDec=decCenter,
                                               boundLength=length,
                                               mjd=51000., bandpassName='i')
            for i in range(5):
                try:
                    dbobj = objcls(verbose=False)
                    result = dbobj.query_columns(column_outputs, obs_metadata=obs_metadata)
                    break
                except:
                    print "network failure. Retrying."

            #testObservationMetadata gives few enough results for one chunk
            try:
                result = result.next()
            except StopIteration:
                raise RuntimeError("No results for %s."%(objname))

            self.assertLess(max(result['raJ2000']), numpy.radians(raCenter+length))
            self.assertGreater(min(result['raJ2000']), numpy.radians(raCenter-length))

            self.assertLess(max(result['decJ2000']), numpy.radians(decCenter+length))
            self.assertGreater(max(result['decJ2000']), numpy.radians(decCenter-length))


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(testCatalogBounds)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
