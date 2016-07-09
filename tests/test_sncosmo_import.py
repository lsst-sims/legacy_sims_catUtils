import unittest
import lsst.utils.tests as utilsTests
import astropy

class sncosmo_import_tests(unittest.TestCase):
    
    def Setup(self):
        pass
    def test_versions(self):
        print(astropy.__version__)
        self.assert(astropy.__version__=='1.2.0')


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(sncosmo_import_tests)
    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)

if __name__ == '__main__':
    run(True)
