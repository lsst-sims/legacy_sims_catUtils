import unittest
import lsst.utils.tests as utilsTests

from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.exampleCatalogDefinitions import ObsStarCatalogBase
#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
import os, inspect

class TestCat(InstanceCatalog):
    catalog_type = 'unit_test_catalog'
    column_outputs = ['raJ2000', 'decJ2000']


class basicAccessTest(unittest.TestCase):
    @unittest.expectedFailure
    def testObjects(self):
        for objname, objcls in CatalogDBObject.registry.iteritems():
            if not objcls.doRunTest or (objcls.testObservationMetaData is None):
                continue
            dbobj = objcls(verbose=False)
            obs_metadata = dbobj.testObservationMetaData
            print "Running tests for", objname
            #Get results all at once
            result = dbobj.query_columns(obs_metadata=obs_metadata)

            #Since there is only one chunck,
            try:
                result = result.next()
            except StopIteration:
                raise RuntimeError("No results for %s defined in %s"%(objname, 
                       inspect.getsourcefile(dbobj.__class__)))
            if objname.startswith('galaxy'):
                TestCat.column_outputs = ['galid', 'raJ2000', 'decJ2000']
            else:
                TestCat.column_outputs = ['raJ2000', 'decJ2000']
            cat = dbobj.getCatalog('unit_test_catalog', obs_metadata)
            if os.path.exists('testCat.out'):
                os.unlink('testCat.out')
            try:
                cat.write_catalog('testCat.out')
            finally:
                if os.path.exists('testCat.out'):
                    os.unlink('testCat.out')

    @unittest.expectedFailure
    def testObsCat(self):
        objname = 'wdstars'
        dbobj = CatalogDBObject.from_objid(objname)
        obs_metadata = dbobj.testObservationMetaData
        # To cover the central ~raft
        obs_metadata.boundLength = 0.4
        opsMetadata = {'Opsim_rotskypos':(0., float),
                       'Unrefracted_RA':(numpy.radians(obs_metadata.unrefractedRA), float),
                       'Unrefracted_Dec':(numpy.radians(obs_metadata.unrefractedDec), float)}
        obs_metadata.phoSimMetaData = opsMetadata
        cat = dbobj.getCatalog('obs_star_cat', obs_metadata)
        if os.path.exists('testCat.out'):
            os.unlink('testCat.out')
        try:
            cat.write_catalog('testCat.out')
        finally:
            if os.path.exists('testCat.out'):
                os.unlink('testCat.out')

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(basicAccessTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
