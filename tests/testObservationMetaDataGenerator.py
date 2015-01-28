import os
import unittest
import numpy
import lsst.utils.tests as utilsTests
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

class ObservationMetaDataGeneratorTest(unittest.TestCase):

    def testQuery(self):
        myGen = ObservationMetaDataGenerator()
        results = myGen.getObservationMetaData(fieldRA=(0.0, numpy.degrees(0.1)))

        print results
        print results['fieldRA']
        for el in results:
            print el['expDate'], el['fiveSigmaDepth']

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(ObservationMetaDataGeneratorTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
