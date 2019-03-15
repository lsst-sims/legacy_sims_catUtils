import unittest
import numpy as np
import time

import lsst.utils.tests


def setup_module(module):
    lsst.utils.tests.init()

class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
