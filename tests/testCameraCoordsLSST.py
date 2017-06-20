"""
This file does not include the boilerplate memory leak tests.  Those tests
kept failing because the instantiations of the LSST Camera were not cleaned
up before the tests were run.
"""

from __future__ import with_statement
import os
import numpy as np
import unittest
import lsst.utils.tests

from lsst.utils import getPackageDir
from lsst.sims.utils import radiansFromArcsec, ObservationMetaData
from lsst.sims.catalogs.db import fileDBObject
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import CameraCoords, CameraCoordsLSST
from lsst.sims.catUtils.mixins import AstrometryStars
from lsst.obs.lsstSim import LsstSimMapper

def setup_module(module):
    lsst.utils.tests.init()


class CameraCoordLSST_testCase(unittest.TestCase):
    """
    This class will test that the CameraCoordsLSST mixin returns
    the same results as the CameraCoords mixin with self.camera=lsst_camera
    """

    def test_different_cameras(self):
        scratch_dir = os.path.join(getPackageDir('sims_catUtils'),'tests',
                                   'scratchSpace')

        db_text_file = os.path.join(scratch_dir, 'cameraCoord_db_text.txt')

        rng = np.random.RandomState(6512)

        pointing_ra = 15.0
        pointing_dec = 13.0

        n_obj = 100
        ra_list = pointing_ra + 2.0*rng.random_sample(n_obj)
        dec_list = pointing_dec + 2.0*rng.random_sample(n_obj)
        px_list = radiansFromArcsec(0.005)*rng.random_sample(n_obj)
        px_list += radiansFromArcsec(0.001)
        mura_list = radiansFromArcsec(0.005)*rng.random_sample(n_obj)
        mudec_list = radiansFromArcsec(0.005)*rng.random_sample(n_obj)
        vrad_list = 100.0*rng.random_sample(n_obj)

        with open(db_text_file, 'w') as out_file:
            for ix, (rdeg, ddeg, rrad, drad, px, mura, mudec, vrad) in \
            enumerate(zip(ra_list, dec_list,
                          np.radians(ra_list), np.radians(dec_list),
                          px_list, mura_list, mudec_list, vrad_list)):

                out_file.write('%d %e %e %e %e %e %e %e %e\n' % (ix, rdeg, ddeg,
                                                                 rrad, drad,
                                                                 px,
                                                                 mura, mudec,
                                                                 vrad))

        dtype = np.dtype([('id', int), ('raDeg', float), ('decDeg', float),
                          ('raJ2000', float), ('decJ2000', float),
                          ('parallax', float),
                          ('properMotionRa', float), ('properMotionDec', float),
                          ('radialVelocity', float)])


        db = fileDBObject(db_text_file, dtype=dtype, idColKey='id')
        db.raColName = 'raDeg'
        db.decColName = 'decDeg'

        if os.path.exists(db_text_file):
            os.unlink(db_text_file)


        class CameraCoordsCatalog(AstrometryStars, CameraCoords,
                                  InstanceCatalog):
            camera = LsstSimMapper().camera
            column_outputs = ['id', 'chipName']


        class CameraCoordsLSSTCatalog(AstrometryStars, CameraCoordsLSST,
                                      InstanceCatalog):
            column_outputs = ['id', 'chipName']

        obs = ObservationMetaData(pointingRA=pointing_ra,
                                  pointingDec=pointing_dec,
                                  boundLength=1.75,
                                  boundType='circle',
                                  rotSkyPos=23.0,
                                  mjd=59580.0)

        control_cat = CameraCoordsCatalog(db, obs_metadata=obs)
        test_cat = CameraCoordsLSSTCatalog(db, obs_metadata=obs)

        control_line_list = []
        none_chips = 0
        for line in control_cat.iter_catalog():
            if line[1] is None:
                none_chips += 1
            control_line_list.append(line)
        self.assertGreater(len(control_line_list), 0)
        self.assertLess(none_chips, len(control_line_list)/2)

        line_ct = 0
        for line in test_cat.iter_catalog():
            line_ct += 1
            self.assertIn(line, control_line_list)
        self.assertEqual(line_ct, len(control_line_list))


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
