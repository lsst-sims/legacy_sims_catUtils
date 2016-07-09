import unittest
import lsst.utils.tests as utilsTests
import astropy
import sncosmo

class sncosmo_import_tests(unittest.TestCase):
    
    def Setup(self):
        import astropy
    def test_versions(self):
        self.assertGreaterEqual(astropy.__version__, '1.1.1')
    def test_versions(self):
        self.assertGreaterEqual(sncosmo.__version__, '1.3.0')


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(sncosmo_import_tests)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
