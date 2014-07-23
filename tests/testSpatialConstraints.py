import unittest
import numpy
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.generation.db import DBObject, ObservationMetaData

#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm

def distanceToPoint(lon1, lat1, lon2, lat2):
    nx1 = numpy.cos(lat1) * numpy.cos(lon1)
    ny1 = numpy.cos(lat1) * numpy.sin(lon1)
    nz1 = numpy.sin(lat1)
    nx2 = numpy.cos(lat2) * numpy.cos(lon2)
    ny2 = numpy.cos(lat2) * numpy.sin(lon2)
    nz2 = numpy.sin(lat2)
    return numpy.rad2deg(2 * numpy.arcsin(
                             numpy.sqrt((nx1 - nx2) * (nx1 - nx2)
                                           + (ny1 - ny2) * (ny1 - ny2)
                                           + (nz1 - nz2) * (nz1 - nz2)
                                        ) / 2.
                                        )
                        )   

vDistanceBetweenPoints = numpy.vectorize(distanceToPoint)

class testCatalogBounds(unittest.TestCase):

    def testCircleBounds(self):
        column_outputs = ['raJ2000', 'decJ2000']
        for objname, objcls in DBObject.registry.iteritems():
            if not objcls.doRunTest \
            or (objcls.testObservationMetaData is None) \
            or (objcls.testObservationMetaData.circ_bounds is None):
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
            self.assertGreater(obs_metadata.circ_bounds['radius'] + 1.e-4, 
                           max(vDistanceBetweenPoints(numpy.radians(obs_metadata.circ_bounds['ra']), 
                                                      numpy.radians(obs_metadata.circ_bounds['dec']),
                                                    result['raJ2000'], result['decJ2000'])))

    def testBoxBounds(self):
        column_outputs = ['raJ2000', 'decJ2000']
        for objname, objcls in DBObject.registry.iteritems():
            if not objcls.doRunTest \
            or (objcls.testObservationMetaData is None) \
            or (objcls.testObservationMetaData.circ_bounds is None):
                continue
            print "Running tests for", objname
            circ_bounds = objcls.testObservationMetaData.circ_bounds 
            obs_metadata = ObservationMetaData(circ_bounds=None, box_bounds=dict(ra_min=circ_bounds['ra'],  
                                               ra_max=circ_bounds['ra'] + 2*circ_bounds['radius'], 
                                               dec_min=circ_bounds['dec'], 
                                               dec_max=circ_bounds['dec'] + 2*circ_bounds['radius']),
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

            self.assertLess(max(result['raJ2000']), numpy.radians(obs_metadata.box_bounds['ra_max']))
            self.assertGreater(min(result['raJ2000']), numpy.radians(obs_metadata.box_bounds['ra_min']))

            self.assertLess(max(result['decJ2000']), numpy.radians(obs_metadata.box_bounds['dec_max']))
            self.assertGreater(max(result['decJ2000']), numpy.radians(obs_metadata.box_bounds['dec_min']))


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
