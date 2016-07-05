import unittest
import lsst.utils.tests as utilsTests

class sncosmo_import_tests(unittest.TestCase):
    
    def Setup(self):
        import astropy
    def test_versions(self):
        print(astropy.__version__)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(sncosmo_import_tests)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
