from builtins import zip
from builtins import range
import numpy as np
import unittest
import lsst.utils.tests

from lsst.sims.utils import ObservationMetaData, _observedFromICRS
from lsst.sims.utils import haversine, arcsecFromRadians
from lsst.sims.catUtils.mixins import PhoSimAstrometryBase


def setup_module(module):
    lsst.utils.tests.init()


class DePrecessionTest(unittest.TestCase):

    def test_de_precession(self):
        """
        test de-precession by de-precessing a list of RA, Dec
        and verifying that the distance between the de-precessed
        points is the same as the distance between the input points.

        Also verify that the observed boresite gets de-precessed correctly
        """
        rng = np.random.RandomState(12)
        n_samples = 5
        pra = 34.0
        pdec = 65.0
        obs = ObservationMetaData(pointingRA=pra,
                                  pointingDec=pdec,
                                  mjd=58324.1)

        raObs, decObs = _observedFromICRS(np.array([np.radians(pra)]),
                                          np.array([np.radians(pdec)]),
                                          obs_metadata=obs, epoch=2000.0,
                                          includeRefraction=False)

        ra_list = []
        dec_list = []
        ra_list.append(raObs[0])
        dec_list.append(decObs[0])
        for rr, dd in zip(rng.random_sample(n_samples)*2.0*np.pi,
                          (rng.random_sample(n_samples)-0.5)*np.pi):

            ra_list.append(rr)
            dec_list.append(dd)

        ra_list = np.array(ra_list)
        dec_list = np.array(dec_list)

        raDecTransformed = PhoSimAstrometryBase()._dePrecess(ra_list, dec_list, obs)
        dd = arcsecFromRadians(haversine(np.radians(pra), np.radians(pdec),
                               raDecTransformed[0][0], raDecTransformed[1][0]))
        self.assertLess(dd, 1.0e-6)
        dd0 = arcsecFromRadians(haversine(raObs[0], decObs[0], np.radians(pra), np.radians(pdec)))
        self.assertLess(dd, dd0)

        for ix in range(n_samples+1):
            for iy in range(n_samples+1):
                if ix != iy:
                    dd1 = arcsecFromRadians(haversine(ra_list[ix], dec_list[ix], ra_list[iy], dec_list[iy]))
                    dd2 = arcsecFromRadians(haversine(raDecTransformed[0][ix], raDecTransformed[1][ix],
                                                      raDecTransformed[0][iy], raDecTransformed[1][iy]))

                    self.assertAlmostEqual(dd1, dd2, delta=6)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
