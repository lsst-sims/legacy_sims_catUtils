import os
import unittest
import lsst.utils.tests as utilsTests

from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.exampleCatalogDefinitions import ObsStarCatalogBase
#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
import os, inspect

class TestCat(InstanceCatalog):
    catalog_type = 'unit_test_catalog'
    column_outputs = ['raJ2000', 'decJ2000']

retry = 5

class basicAccessTest(unittest.TestCase):
    @unittest.expectedFailure
    def testObjects(self):
        for objname, objcls in CatalogDBObject.registry.iteritems():
            if not objcls.doRunTest or (objcls.testObservationMetaData is None):
                continue
            for i in range(retry):   
                try:
                    dbobj = objcls(verbose=False)
                    break
                except:
                    continue
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

    def testObsCat(self):
        objname = 'wdstars'
        for i in range(retry):   
            try:
                dbobj = CatalogDBObject.from_objid(objname)
                break
            except:
                continue
        obs_metadata = dbobj.testObservationMetaData
        # To cover the central ~raft
        obs_metadata.boundLength = 0.4
        opsMetadata = {'Opsim_rotskypos':(0., float),
                       'Unrefracted_RA':(obs_metadata.unrefractedRA, float),
                       'Unrefracted_Dec':(obs_metadata.unrefractedDec, float)}
        obs_metadata.phoSimMetadata = opsMetadata
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
