"""
This unit test essentially copies the unit test of the same name that
exists in sims_utils.  However, in this iteration, the test to verify
that a mapped filename exists is never skipped.  The idea being that we need
to make sure that we have loaded the correct version of sims_sed_library.
"""

import os
import unittest
import lsst.utils.tests as utilsTests

from lsst.utils import getPackageDir
from lsst.sims.utils import defaultSpecMap

class FileMapTest(unittest.TestCase):

    def setUp(self):
        self.root_dir = getPackageDir('sims_sed_library')


    def verifyFile(self, file_name, dir_name):
        """
        Verify that specMape[file_name] results in a file in dir dir_name
        """
        test_name = defaultSpecMap[file_name]
        control_name = os.path.join(dir_name, file_name+'.gz')
        msg = '%s should map to %s; it actually maps to %s' % (file_name, control_name, test_name)
        self.assertEqual(test_name, control_name, msg=msg)

        add_space = file_name+' '
        self.assertNotEqual(add_space, file_name)
        test_name = defaultSpecMap[add_space]
        msg = '%s should map to %s; it actually maps to %s' % (add_space, control_name, test_name)
        self.assertEqual(test_name, control_name, msg=msg)

        add_space = ' '+file_name
        self.assertNotEqual(add_space, file_name)
        test_name = defaultSpecMap[add_space]
        msg = '%s should map to %s; it actually maps to %s' % (add_space, control_name, test_name)
        self.assertEqual(test_name, control_name, msg=msg)

        add_gz = file_name+'.gz'
        self.assertNotEqual(add_gz, file_name)
        test_name = defaultSpecMap[add_gz]
        msg = '%s should map to %s; it actually maps to %s' % (add_gz, control_name, test_name)
        self.assertEqual(test_name, control_name, msg=msg)

        # 12 January 2016
        # Currently, this test must be commented out because
        # Jenkins does not re-run eupspkg.cfg.sh every time
        # it tests a branch, which means that it won't get the
        # new library.  Once we move sims_sed_library to a git LFS
        # repo, this test should be performed.
        #
        #full_path = os.path.join(self.root_dir, test_name)
        #msg = '%s does not exist; it should' % full_path
        #self.assertTrue(os.path.exists(full_path), msg=msg)


    def testMLT(self):
        """
        Test that defaultSpecMap correctly locates MLT dwarf spectra
        """
        self.verifyFile('lte004-3.5-0.0a+0.0.BT-Settl.spec', 'starSED/mlt')


    def test_m_spec(self):
        """
        Test that defaultSpecMap correctly finds old MLT dwarf spectra
        that begin with 'm'
        """
        self.verifyFile('m5.1Full.dat', 'starSED/old_mlt')


    def test_l4_spec(self):
        """
        Test that defaultSpecMap correctly finds l4Full.dat
        """
        self.verifyFile('l4Full.dat', 'starSED/old_mlt')


    def test_L_spec(self):
        """
        Test that defaultSpecMap correctly find the L#_# spectra
        """
        self.verifyFile('L2_0Full.dat', 'starSED/old_mlt')


    def test_burrows_spec(self):
        """
        Test that defaultSpecMap correctly find the burrows spectra
        """
        self.verifyFile('burrows+2006c91.21_T1400_g5.5_cf_0.3x', 'starSED/old_mlt')


    def testBergeron(self):
        """
        Test that defaultSpecMap correctly locates the bergeron spectra
        """
        self.verifyFile('bergeron_4750_85.dat_4900', 'starSED/wDs')


    def testKurucz(self):
        """
        Test that defaultSpecMap correctly locates the kurucz spectra
        """
        self.verifyFile('km30_5000.fits_g10_5040', 'starSED/kurucz')
        self.verifyFile('kp10_9000.fits_g40_9100', 'starSED/kurucz')


    def testGalaxy(self):
        """
        Test that defaultSpecMap correctly locates the galaxy SEDs
        """
        self.verifyFile('Const.79E06.002z.spec', 'galaxySED')
        self.verifyFile('Inst.79E06.02Z.spec', 'galaxySED')
        self.verifyFile('Exp.40E08.02Z.spec', 'galaxySED')
        self.verifyFile('Burst.40E08.002z.spec', 'galaxySED')

def suite():
    """Returns a suite containing all the test cases in this module."""
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(FileMapTest)

    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
