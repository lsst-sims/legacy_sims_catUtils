from __future__ import print_function
from builtins import zip
from builtins import next
import os
import inspect
import numpy as np
import sys
import traceback
import unittest
import tempfile
import shutil
import lsst.utils.tests

from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.utils import failedOnFatboy
# The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm


def setup_module(module):
    lsst.utils.tests.init()


class DummyCat(InstanceCatalog):
    catalog_type = __file__ + 'unit_test_catalog'
    column_outputs = ['raJ2000', 'decJ2000']


class basicAccessTest(unittest.TestCase):

    longMessage = True

    def testObjects(self):
        catDir = tempfile.mkdtemp('basicAccessTest_testObjects')
        if not os.path.exists(catDir):
            os.mkdir(catDir)
        catName = tempfile.mktemp(prefix='basicAccessTest_testObjects',
                                  dir=catDir, suffix='.txt')
        ct_connected = 0
        ct_failed_connection = 0
        list_of_failures = []

        for objname, objcls in CatalogDBObject.registry.items():
            if not objcls.doRunTest or (objcls.testObservationMetaData is None):
                continue

            print("Running tests for", objname)
            try:
                dbobj = objcls(verbose=False)
            except:
                trace = traceback.extract_tb(sys.exc_info()[2], limit=20)
                msg = sys.exc_info()[1].args[0]
                if 'Failed to connect' in str(msg) or failedOnFatboy(trace):

                    # if the exception was due to a failed connection
                    # to fatboy, ignore it

                    ct_failed_connection += 1
                    list_of_failures.append(objname)
                    continue
                else:
                    raise

            obs_metadata = dbobj.testObservationMetaData

            # Get results all at once
            try:
                result = dbobj.query_columns(obs_metadata=obs_metadata)
            except:

                # This is because the solar system object 'tables'
                # don't actually connect to tables on fatboy; they just
                # call methods stored on fatboy.  Therefore, the connection
                # failure will not be noticed until this part of the test

                ct_failed_connection += 1
                list_of_failures.append(objname)
                msg = sys.exc_info()[1].args[0]
                if 'DB-Lib error' in msg:
                    continue
                else:
                    raise

            ct_connected += 1

            # Since there is only one chunk,
            try:
                result = next(result)
            except StopIteration:
                raise RuntimeError("No results for %s defined in %s"%(objname,
                                   inspect.getsourcefile(dbobj.__class__)))
            if objname.startswith('galaxy'):
                DummyCat.column_outputs = ['galid', 'raJ2000', 'decJ2000']
            else:
                DummyCat.column_outputs = ['raJ2000', 'decJ2000']
            cat = dbobj.getCatalog(__file__+'unit_test_catalog', obs_metadata)
            if os.path.exists(catName):
                os.unlink(catName)
            try:
                cat.write_catalog(catName)
                dtypeList = [(name, np.float) for name in cat._column_outputs]
                testData = np.genfromtxt(catName, delimiter = ', ',
                                         dtype=np.dtype(dtypeList))
                self.assertGreater(len(testData), 0)
            finally:
                if os.path.exists(catName):
                    os.unlink(catName)

        if os.path.exists(catDir):
            shutil.rmtree(catDir, ignore_errors=True)

        self.assertEqual(len(list_of_failures), ct_failed_connection)

        print('\n================')
        print('Do not worry about this message')
        print('sometimes, connections to the UW database fail.')
        print('It is expected.')
        print('This is just a tally so that you know how often that happened.')
        print('successful connections: ', ct_connected)
        print('failed connections: ', ct_failed_connection)
        if len(list_of_failures) > 0:
            print('objects that failed to connect: ', list_of_failures)

    def testObsCat(self):
        objname = 'wdstars'
        catDir = tempfile.mkdtemp('basicAccessTest_testObsCat')
        if not os.path.exists(catDir):
            os.mkdir(catDir)
        catName = tempfile.mktemp(prefix='basicAccessTest_testObsCat',
                                  dir=catDir, suffix='.txt')

        try:
            dbobj = CatalogDBObject.from_objid(objname)
            obs_metadata = dbobj.testObservationMetaData
            # To cover the central ~raft
            obs_metadata.boundLength = 0.4
            obs_metadata.rotSkyPos = 0.0
            cat = dbobj.getCatalog('obs_star_cat', obs_metadata)
            if os.path.exists(catName):
                os.unlink(catName)
            try:
                cat.write_catalog(catName)
                dtypeList = [(name, np.float) for name in cat._column_outputs]
                testData = np.genfromtxt(catName, delimiter = ', ',
                                         dtype=np.dtype(dtypeList))
                self.assertGreater(len(testData), 0)
            finally:
                if os.path.exists(catName):
                    os.unlink(catName)
                if os.path.exists(catDir):
                    shutil.rmtree(catDir, ignore_errors=True)

            print('\ntestObsCat successfully connected to fatboy')

        except:
            trace = traceback.extract_tb(sys.exc_info()[2], limit=20)
            msg = sys.exc_info()[1].args[0]
            if 'Failed to connect' in str(msg) or failedOnFatboy(trace):

                # if the exception was because of a failed connection
                # to fatboy, ignore it.

                print('\ntestObsCat failed to connect to fatboy')
                print('Sometimes that happens.  Do not worry.')

                if os.path.exists(catDir):
                    shutil.rmtree(catDir, ignore_errors=True)

                pass
            else:
                raise

    def test_limit(self):
        """
        Test that the limit kwarg in query_columns behaves correctly.

        Will test on one star table and one galaxy table.
        """
        list_of_failures = []
        for objcls, clsname in zip((bcm.StarObj, bcm.GalaxyObj), ('StarObj', 'GalaxyObj')):
            msg = "failed the limit test\noffending class is %s" % clsname
            try:
                dbobj = objcls(verbose=False)
            except:
                trace = traceback.extract_tb(sys.exc_info()[2], limit=20)
                msg = sys.exc_info()[1].args[0]
                if 'Failed to connect' in str(msg) or failedOnFatboy(trace):

                    # if the exception was due to a failed connection
                    # to fatboy, ignore it
                    list_of_failures.append(clsname)
                    continue
                else:
                    raise

            obs_metadata = dbobj.testObservationMetaData

            results = dbobj.query_columns(obs_metadata=obs_metadata)

            ct_res = 0
            for chunk in results:
                for line in chunk:
                    ct_res += 1

            self.assertGreater(ct_res, 10, msg=msg)

            limited_results = dbobj.query_columns(obs_metadata=obs_metadata, limit=10)

            ct_limit = 0
            for chunk in limited_results:
                for line in chunk:
                    ct_limit += 1

            self.assertEqual(ct_limit, 10, msg=msg)

        if len(list_of_failures) > 0:
            print("\nList of DBObjects that could not connect to fatboy " \
                  "for the test on the limit kwarg")
            for nn in list_of_failures:
                print(nn)

    def test_constraint(self):
        """
        Test that passing a constraint into query_columns works (i.e. if I only want
        to select galaxies whose varParamStr is not NULL).
        """
        list_of_failures = []
        constraint = "varParamStr IS NOT NULL"
        for objcls, clsname in zip((bcm.GalaxyObj, bcm.GalaxyTileObj), ('GalaxyObj', 'GalaxyTileObj')):
            msg = "failed the constraint test\noffending class is %s" % clsname
            try:
                dbobj = objcls(verbose=False)
            except:
                trace = traceback.extract_tb(sys.exc_info()[2], limit=20)
                msg = sys.exc_info()[1].args[0]
                if 'Failed to connect' in str(msg) or failedOnFatboy(trace):

                    # if the exception was due to a failed connection
                    # to fatboy, ignore it
                    list_of_failures.append(clsname)
                    continue
                else:
                    raise

            obs_metadata = dbobj.testObservationMetaData

            # query witout a constraint on varParamStr
            results = dbobj.query_columns(colnames=['raJ2000', 'decJ2000', 'varParamStr'],
                                          obs_metadata=obs_metadata)

            # count total number of rows (ct_res) and number of rows with a null
            # varParamStr (ct_no_varparamstr).  Note that varParamStr will be the
            # index=3 entry in result rows because the id of the galaxy gets
            # automatically added to query results.
            ct_res = 0
            ct_no_varparamstr = 0
            for chunk in results:
                for line in chunk:
                    if line[3] == 'None':
                        ct_no_varparamstr += 1
                    ct_res += 1

            # run the same query, but demanding that varParamStr is not NULL
            constrained_results = dbobj.query_columns(colnames=['raJ2000', 'decJ2000', 'varParamStr'],
                                                      obs_metadata=obs_metadata,
                                                      constraint=constraint)

            # count the number of rows with non-NULL varParamStr
            ct_con = 0
            for chunk in constrained_results:
                for line in chunk:
                    ct_con += 1
                    self.assertNotEqual(line[3], 'None')

            # check that the number of non-NULL varParamStr and NULL varParamStr rows
            # compare the way that they should
            self.assertGreater(ct_res, ct_con)
            self.assertGreater(ct_no_varparamstr, 0)
            self.assertEqual(ct_res-ct_con, ct_no_varparamstr)

        if len(list_of_failures) > 0:
            print("\nList of DBObjects that could not connect to fatboy " \
                  "for the test on the constraint kwarg")
            for nn in list_of_failures:
                print(nn)

    def test_limit_and_constraint(self):
        """
        Test that limit and constraint work together
        """
        list_of_failures = []
        constraint = "varParamStr IS NOT NULL"
        could_connect = True
        try:
            dbobj = bcm.GalaxyObj(verbose=False)
        except:
            trace = traceback.extract_tb(sys.exc_info()[2], limit=20)
            msg = sys.exc_info()[1].args[0]
            if 'Failed to connect' in str(msg) or failedOnFatboy(trace):

                # if the exception was due to a failed connection
                # to fatboy, ignore it
                list_of_failures.append('GalaxyObj')
                could_connect = False
            else:
                raise

        if could_connect:
            obs_metadata = dbobj.testObservationMetaData

            # query with a constraint on varParamStr but no limit
            results_no_limit = dbobj.query_columns(colnames=['raJ2000', 'decJ2000', 'varParamStr'],
                                                   obs_metadata=obs_metadata,
                                                   constraint=constraint)

            ct_res = 0
            for chunk in results_no_limit:
                for line in chunk:
                    self.assertNotEqual(line[3], 'None')
                    ct_res += 1

            self.assertGreater(ct_res, 1)

            # run the same query, but limiting the results
            limited_results = dbobj.query_columns(colnames=['raJ2000', 'decJ2000', 'varParamStr'],
                                                  obs_metadata=obs_metadata,
                                                  constraint=constraint,
                                                  limit=ct_res-1)
            ct_lim = 0
            for chunk in limited_results:
                for line in chunk:
                    ct_lim += 1
                    self.assertNotEqual(line[3], 'None')

                self.assertEqual(ct_lim, ct_res-1)

        if len(list_of_failures) > 0:
            print("\nList of DBObjects that could not connect to fatboy " \
                  "for the test on the constraint and limit kwargs")
            for nn in list_of_failures:
                print(nn)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
