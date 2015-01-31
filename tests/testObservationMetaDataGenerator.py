from __future__ import with_statement
import os
import unittest
import numpy
import lsst.utils.tests as utilsTests
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catalogs.generation.db import CatalogDBObject, CircleBounds, BoxBounds
from lsst.sims.catUtils.baseCatalogModels import GalaxyBulgeObj, GalaxyDiskObj, GalaxyAgnObj, StarObj
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB

#28 January 2015
#SearchReversion and testGalaxyBulge are duplicated from testPhoSimCatalogs.py
#There is already a pull request outstanding that will move
#these classes to a testUtils.py file in sims_catUtils/../utils/
#Once that pull request is merged, I will clean this file up
#and remove the duplicated code.
class SearchReversion(CatalogDBObject):
    """
    This is a mixin which is used below to force the galaxy CatalogDBObjects created for
    this unittest to use the methods defined in the CatalogDBObject class.  This is because
    we are using classes which inherit from GalaxyTileObj but do not actually want to use
    the tiled query routines.

    We also use this mixin for our stellar database object.  This is because StarObj
    implements a query search based on htmid, which the test database for this unit
    test will not have.
    """

    def _get_column_query(self, *args, **kwargs):
        return CatalogDBObject._get_column_query(self,*args, **kwargs)

    def _final_pass(self, *args, **kwargs):
        return CatalogDBObject._final_pass(self,*args, **kwargs)

    def query_columns(self, *args, **kwargs):
        return CatalogDBObject.query_columns(self, *args, **kwargs)

class testGalaxyBulge(SearchReversion, GalaxyBulgeObj):
    """
    A class for storing galaxy bulges
    """
    objid = 'phoSimTestBulges'
    objectTypeId = 88

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = GalaxyBulgeObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class ObservationMetaDataGeneratorTest(unittest.TestCase):

    def testExceptions(self):
        """
        Make sure that RuntimeErrors get raised when they should
        """
        gen = ObservationMetaDataGenerator()
        self.assertRaises(RuntimeError, gen.getObservationMetaData)
        self.assertRaises(RuntimeError, gen.getObservationMetaData,fieldRA=(1.0, 2.0, 3.0))

    def testQueryOnRanges(self):
        """
        Test that ObservationMetaData objects returned by queries of the form
        min < value < max
        are, in fact, within that range.

        Test when querying on both a single and two columns.
        """
        gen = ObservationMetaDataGenerator()

        #An list containing the bounds of our queries.
        #The order of the tuples must correspond to the order of
        #self.columnMapping in ObservationMetaDataGenerator.
        #This was generated with a separate script which printed
        #the median and maximum values of all of the quantities
        #in our test opsim database
        bounds = [
        ('obsHistID',(5973, 11080)),
        ('expDate',(1220779, 1831593)),
        ('fieldRA',(numpy.degrees(1.370916), numpy.degrees(1.5348635))),
        ('fieldDec',(numpy.degrees(-0.456238), numpy.degrees(0.0597905))),
        ('moonRA',(numpy.degrees(2.914132), numpy.degrees(4.5716525))),
        ('moonDec',(numpy.degrees(0.06305), numpy.degrees(0.2216745))),
        ('rotSkyPos',(numpy.degrees(3.116656), numpy.degrees(4.6974265))),
        ('telescopeFilter',('i','i')),
        ('rawSeeing',(0.728562, 1.040495)),
        ('sunAlt',(numpy.degrees(-0.522905), numpy.degrees(-0.366073))),
        ('moonAlt',(numpy.degrees(0.099096), numpy.degrees(0.5495415))),
        ('dist2Moon',(numpy.degrees(1.570307), numpy.degrees(2.347868))),
        ('moonPhase',(52.2325, 76.0149785)),
        ('expMJD',(49367.129396, 49374.1990025)),
        ('altitude',(numpy.degrees(0.781015), numpy.degrees(1.1433785))),
        ('azimuth',(numpy.degrees(3.470077), numpy.degrees(4.8765995))),
        ('visitExpTime',(30.0,30.0)),
        ('airmass',(1.420459, 2.0048075)),
        ('skyBrightness',(19.017605, 20.512553)),
        ('m5',(22.815249, 24.0047695))]

        #test querying on a single column
        for (ii,line) in enumerate(bounds):
            tag = line[0]
            if tag != 'telescopeFilter' and tag != 'visitExpTime':
                args = {}
                args[tag] = line[1]
                results = gen.getObservationMetaData(**args)

                name = gen.columnMapping[ii][2]
                if name is not None:
                    if gen.columnMapping[ii][4] is not None:
                        xmin = gen.columnMapping[ii][4](line[1][0])
                        xmax = gen.columnMapping[ii][4](line[1][1])
                    else:
                        xmin = line[1][0]
                        xmax = line[1][1]
                    ct = 0
                    for obs_metadata in results:
                        ct += 1
                        self.assertTrue(obs_metadata.phoSimMetadata[name][0]<xmax)
                        self.assertTrue(obs_metadata.phoSimMetadata[name][0]>xmin)
                        if tag == 'm5':
                            self.assertTrue(obs_metadata.m5('i')<xmax)
                            self.assertTrue(obs_metadata.m5('i')>xmin)

                    #make sure that we did not accidentally choose values such that
                    #no ObservationMetaData were ever returned
                    self.assertTrue(ct>0)

        #test querying on two columns at once
        ct = 0
        for ii in range(len(bounds)):
            tag1 = bounds[ii][0]
            if tag1 != 'telescopeFilter' and tag1 != 'visitExpTime':
                name1 = gen.columnMapping[ii][2]
                if gen.columnMapping[ii][4] is not None:
                    xmin = gen.columnMapping[ii][4](bounds[ii][1][0])
                    xmax = gen.columnMapping[ii][4](bounds[ii][1][1])
                else:
                    xmin = bounds[ii][1][0]
                    xmax = bounds[ii][1][1]
                for jj in range(ii+1, len(bounds)):
                    tag2 = bounds[jj][0]
                    if tag2 != 'telescopeFilter' and tag2 != 'visitExpTime':
                        name2 = gen.columnMapping[jj][2]
                        if gen.columnMapping[jj][4] is not None:
                            ymin = gen.columnMapping[jj][4](bounds[jj][1][0])
                            ymax = gen.columnMapping[jj][4](bounds[jj][1][1])
                        else:
                            ymin = bounds[jj][1][0]
                            ymax = bounds[jj][1][1]
                        args = {}
                        args[tag1] = bounds[ii][1]
                        args[tag2] = bounds[jj][1]
                        results = gen.getObservationMetaData(**args)
                        if name1 is not None or name2 is not None:
                            for obs_metadata in results:
                                ct += 1
                                if name1 is not None:
                                    self.assertTrue(obs_metadata.phoSimMetadata[name1][0]>xmin)
                                    self.assertTrue(obs_metadata.phoSimMetadata[name1][0]<xmax)
                                if name2 is not None:
                                    self.assertTrue(obs_metadata.phoSimMetadata[name2][0]>ymin)
                                    self.assertTrue(obs_metadata.phoSimMetadata[name2][0]<ymax)

        #Make sure that we didn't choose values such that no ObservationMetaData were
        #ever returned
        self.assertTrue(ct>0)

    def testQueryExactValues(self):
        """
        Test that ObservationMetaData returned by a query demanding an exact value do,
        in fact, adhere to that requirement.
        """
        gen = ObservationMetaDataGenerator()

        bounds = [
        ('obsHistID',5973),
        ('expDate',1220779),
        ('fieldRA',numpy.degrees(1.370916)),
        ('fieldDec',numpy.degrees(-0.456238)),
        ('moonRA',numpy.degrees(2.914132)),
        ('moonDec',numpy.degrees(0.06305)),
        ('rotSkyPos',numpy.degrees(3.116656)),
        ('telescopeFilter','i'),
        ('rawSeeing',0.728562),
        ('sunAlt',numpy.degrees(-0.522905)),
        ('moonAlt',numpy.degrees(0.099096)),
        ('dist2Moon',numpy.degrees(1.570307)),
        ('moonPhase',52.2325),
        ('expMJD',49367.129396),
        ('altitude',numpy.degrees(0.781015)),
        ('azimuth',numpy.degrees(3.470077)),
        ('visitExpTime',30.0),
        ('airmass',1.420459),
        ('skyBrightness',19.017605),
        ('m5',22.815249)]

        for ii in range(len(bounds)):
            tag = bounds[ii][0]
            if tag != 'telescopeFilter' and tag != 'visitExpTime':
                name = gen.columnMapping[ii][2]
                args = {}
                args[tag] = bounds[ii][1]
                results = gen.getObservationMetaData(**args)

                if gen.columnMapping[ii][4] is not None:
                    value = gen.columnMapping[ii][4](bounds[ii][1])
                else:
                    value = bounds[ii][1]

                if name is not None:
                    ct = 0
                    for obs_metadata in results:
                        self.assertAlmostEqual(value, obs_metadata.phoSimMetadata[name][0],10)
                        ct += 1

                    #Make sure that we did not choose a value which returns zero ObservationMetaData
                    self.assertTrue(ct>0)

    def testQueryLimit(self):
        """
        Test that, when we specify a limit on the number of ObservationMetaData we want returned,
        that limit is respected
        """
        gen = ObservationMetaDataGenerator()
        results = gen.getObservationMetaData(fieldRA=(numpy.degrees(1.370916), numpy.degrees(1.5348635)),
                                             limit=20)
        self.assertEqual(len(results),20)

    def testQueryOnFilter(self):
        """
        Test that queries on the filter work.
        """
        gen = ObservationMetaDataGenerator()
        results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916), telescopeFilter='i')
        ct = 0
        for obs_metadata in results:
            self.assertAlmostEqual(obs_metadata.phoSimMetadata['Unrefracted_RA'][0],1.370916)
            self.assertEqual(obs_metadata.phoSimMetadata['Opsim_filter'][0],'i')
            ct += 1

        #Make sure that more than zero ObservationMetaData were returned
        self.assertTrue(ct>0)

    def testObsMetaDataBounds(self):
        """
        Make sure that the bound specifications (i.e. a circle or a box on the sky) are correctly
        passed through to the resulting ObservationMetaData
        """

        gen = ObservationMetaDataGenerator()

        #Test a cirlce with a specified radius
        results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916), telescopeFilter='i', boundLength=0.9)
        ct = 0
        for obs_metadata in results:
            self.assertTrue(isinstance(obs_metadata.bounds,CircleBounds))
            self.assertAlmostEqual(obs_metadata.bounds.radiusdeg,0.9,10)
            self.assertAlmostEqual(obs_metadata.bounds.RA,obs_metadata.phoSimMetadata['Unrefracted_RA'][0],10)
            self.assertAlmostEqual(obs_metadata.bounds.DEC,obs_metadata.phoSimMetadata['Unrefracted_Dec'][0],10)
            ct += 1

        #Make sure that some ObservationMetaData were tested
        self.assertTrue(ct>0)

        #test a square
        results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916), telescopeFilter='i',
                                             boundType='box', boundLength=1.2)
        ct = 0
        for obs_metadata in results:
            RAdeg = numpy.degrees(obs_metadata.phoSimMetadata['Unrefracted_RA'][0])
            DECdeg = numpy.degrees(obs_metadata.phoSimMetadata['Unrefracted_Dec'][0])
            self.assertTrue(isinstance(obs_metadata.bounds,BoxBounds))

            self.assertAlmostEqual(obs_metadata.bounds.RAminDeg,RAdeg-1.2,10)
            self.assertAlmostEqual(obs_metadata.bounds.RAmaxDeg,RAdeg+1.2,10)
            self.assertAlmostEqual(obs_metadata.bounds.DECminDeg,DECdeg-1.2,10)
            self.assertAlmostEqual(obs_metadata.bounds.DECmaxDeg,DECdeg+1.2,10)

            self.assertAlmostEqual(obs_metadata.bounds.RA,obs_metadata.phoSimMetadata['Unrefracted_RA'][0],10)
            self.assertAlmostEqual(obs_metadata.bounds.DEC,obs_metadata.phoSimMetadata['Unrefracted_Dec'][0],10)

            ct += 1

        #Make sure that some ObservationMetaData were tested
        self.assertTrue(ct>0)

        #test a rectangle
        results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916), telescopeFilter='i',
                                             boundType='box', boundLength=(1.2,0.6))
        ct = 0
        for obs_metadata in results:
            RAdeg = numpy.degrees(obs_metadata.phoSimMetadata['Unrefracted_RA'][0])
            DECdeg = numpy.degrees(obs_metadata.phoSimMetadata['Unrefracted_Dec'][0])
            self.assertTrue(isinstance(obs_metadata.bounds,BoxBounds))

            self.assertAlmostEqual(obs_metadata.bounds.RAminDeg,RAdeg-1.2,10)
            self.assertAlmostEqual(obs_metadata.bounds.RAmaxDeg,RAdeg+1.2,10)
            self.assertAlmostEqual(obs_metadata.bounds.DECminDeg,DECdeg-0.6,10)
            self.assertAlmostEqual(obs_metadata.bounds.DECmaxDeg,DECdeg+0.6,10)

            self.assertAlmostEqual(obs_metadata.bounds.RA,obs_metadata.phoSimMetadata['Unrefracted_RA'][0],10)
            self.assertAlmostEqual(obs_metadata.bounds.DEC,obs_metadata.phoSimMetadata['Unrefracted_Dec'][0],10)

            ct += 1

        #Make sure that some ObservationMetaData were tested
        self.assertTrue(ct>0)

    def testCreationOfPhoSimCatalog(self):
        """
        Make sure that we can create PhoSim input catalogs using the returned ObservationMetaData.
        This test will just make sure that all of the expected header entries are there.
        """

        dbName = 'obsMetaDataGeneratorTest.db'
        catName = 'testPhoSimFromObsMetaDataGenerator.txt'
        if os.path.exists(dbName):
            os.unlink(dbName)
        junk_obs_metadata = makePhoSimTestDB(filename=dbName)
        bulgeDB = testGalaxyBulge(address='sqlite:///'+dbName)

        gen = ObservationMetaDataGenerator()
        results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916),telescopeFilter='i')
        testCat = PhoSimCatalogSersic2D(bulgeDB, obs_metadata=results[0])
        testCat.write_catalog(catName)

        filterTranslation=['u','g','r','i','z','y']

        with open(catName) as inputFile:
            lines = inputFile.readlines()
            ix = 0
            for control in gen.columnMapping:
                if control[0] != 'm5' and control[0]!='skyBrightness':
                    words = lines[ix].split()
                    self.assertEqual(control[2], words[0])

                    #skip telescope filter since it doesn't have a mapping in ObservationMetaDataGenerator
                    if control[0] != 'telescopeFilter':
                        if control[4] is not None:
                            value = control[4](float(words[1]))
                        else:
                            value = float(words[1])

                        self.assertAlmostEqual(value, results[0].phoSimMetadata[control[2]][0], 5)
                    else:
                        self.assertEqual(filterTranslation[int(words[1])],results[0].phoSimMetadata[control[2]][0])
                        
                    ix += 1

        if os.path.exists(catName):
            os.unlink(catName)

        if os.path.exists(dbName):
            os.unlink(dbName)



def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(ObservationMetaDataGeneratorTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
