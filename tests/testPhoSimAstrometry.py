import unittest
import os
import numpy as np
import tempfile
import palpy

import lsst.utils.tests
from lsst.utils import getPackageDir

from lsst.sims.utils import _angularSeparation
from lsst.sims.utils import angularSeparation
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.utils import makePhoSimTestDB
from lsst.sims.catUtils.utils import testStarsDBObj, testGalaxyDiskDBObj
from lsst.sims.catUtils.mixins import PhoSimAstrometryBase
from lsst.sims.catUtils.mixins import PhoSimAstrometryStars
from lsst.sims.catUtils.mixins import PhoSimAstrometryGalaxies

from lsst.sims.utils import _observedFromICRS
from lsst.sims.utils import _observedFromAppGeo


from lsst.sims.utils import observedFromICRS
from lsst.sims.utils import observedFromAppGeo

from lsst.sims.utils import Site, ObservationMetaData
from lsst.sims.utils import distanceToSun
from lsst.sims.utils import raDecFromAltAz
from lsst.sims.utils import pupilCoordsFromRaDec
from lsst.sims.coordUtils import focalPlaneCoordsFromPupilCoords
from lsst.sims.coordUtils import lsst_camera

from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.utils.CodeUtilities import _validate_inputs

ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()


def _naivePupilCoordsFromObserved(ra_obs, dec_obs, ra0, dec0, rotSkyPos):
    """
    Convert Observed RA, Dec into pupil coordinates

    Parameters
    ----------
    ra_obs is the observed RA in radians

    dec_obs is the observed Dec in radians

    obs_metadata is an ObservationMetaData characterizing the telescope location and pointing

    epoch is the epoch of the mean RA and Dec in julian years (default=2000.0)

    includeRefraction is a boolean controlling the application of refraction.

    Returns
    --------
    A numpy array whose first row is the x coordinate on the pupil in
    radians and whose second row is the y coordinate in radians
    """

    are_arrays = _validate_inputs([ra_obs, dec_obs], ['ra_obs', 'dec_obs'],
                                  "pupilCoordsFromObserved")

    theta = -1.0*rotSkyPos

    ra_pointing = ra0
    dec_pointing = dec0

    # palpy.ds2tp performs the gnomonic projection on ra_in and dec_in
    # with a tangent point at (pointingRA, pointingDec)
    #
    if not are_arrays:
        try:
            x, y = palpy.ds2tp(ra_obs, dec_obs, ra_pointing, dec_pointing)
        except:
            x = np.NaN
            y = np.NaN
    else:
        try:
            x, y = palpy.ds2tpVector(ra_obs, dec_obs, ra_pointing, dec_pointing)
        except:
            # apparently, one of your ra/dec values was improper; we will have to do this
            # element-wise, putting NaN in the place of the bad values
            x = []
            y = []
            for rr, dd in zip(ra_obs, dec_obs):
                try:
                    xx, yy = palpy.ds2tp(rr, dd, ra_pointing, dec_pointing)
                except:
                    xx = np.NaN
                    yy = np.NaN
                x.append(xx)
                y.append(yy)
            x = np.array(x)
            y = np.array(y)

    # rotate the result by rotskypos (rotskypos being "the angle of the sky relative to
    # camera coordinates" according to phoSim documentation) to account for
    # the rotation of the focal plane about the telescope pointing

    x_out = x*np.cos(theta) - y*np.sin(theta)
    y_out = x*np.sin(theta) + y*np.cos(theta)

    return np.array([x_out, y_out])


class StarTestCatalog(PhoSimAstrometryStars, InstanceCatalog):
    column_outputs = ['raICRS', 'decICRS', 'raPhoSim', 'decPhoSim',
                      'raJ2000', 'decJ2000',
                      'properMotionRa','properMotionDec', 'parallax',
                      'radialVelocity']

    transformations = {'raICRS': np.degrees, 'decICRS': np.degrees,
                       'raPhoSim': np.degrees, 'decPhoSim': np.degrees}

    default_formats = {'f': '%.12g'}

    delimiter = ' '


class GalaxyTestCatalog(PhoSimAstrometryGalaxies, InstanceCatalog):
    column_outputs = ['raICRS', 'decICRS', 'raPhoSim', 'decPhoSim']
    transformations = {'raICRS': np.degrees, 'decICRS': np.degrees,
                       'raPhoSim': np.degrees, 'decPhoSim': np.degrees}

    override_formats = {'raICRS': '%.12g', 'decICRS': '%.12g',
                        'raPhoSim': '%.12g', 'decPhoSim': '%.12g'}

    delimiter = ' '


class PhoSimAstrometryTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.db_name = tempfile.mktemp(dir=ROOT, prefix='PhoSimAstDB', suffix='.db')
        cls.obs = makePhoSimTestDB(filename=cls.db_name,
                                   size=1000)

    @classmethod
    def tearDownClass(cls):
        if hasattr(lsst_camera, '_lsst_camera'):
            del lsst_camera._lsst_camera

        sims_clean_up()
        if os.path.exists(cls.db_name):
            os.unlink(cls.db_name)

    def test_stellar_astrometry_radians(self):
        """
        Test that we can go from raPhoSim, decPhoSim to ICRS coordinates
        in the case of stars (in radians)
        """
        db = testStarsDBObj(driver='sqlite', database=self.db_name)
        cat = StarTestCatalog(db, obs_metadata=self.obs)
        with lsst.utils.tests.getTempFilePath('.txt') as cat_name:
            cat.write_catalog(cat_name)
            dtype = np.dtype([('raICRS', float), ('decICRS', float),
                             ('raPhoSim', float), ('decPhoSim', float),
                             ('raJ2000', float), ('decJ2000', float),
                             ('pmRA', float), ('pmDec', float),
                             ('parallax', float), ('vRad', float)])
            data = np.genfromtxt(cat_name, dtype=dtype)
        self.assertGreater(len(data), 100)
        ra_pho_rad = np.radians(data['raPhoSim'])
        dec_pho_rad = np.radians(data['decPhoSim'])

        # verify that, when transforming back to ICRS, we are within
        # 10^-3 arcsec
        ra_icrs, dec_icrs = PhoSimAstrometryBase._icrsFromPhoSim(ra_pho_rad,
                                                                 dec_pho_rad,
                                                                 self.obs)
        dist = _angularSeparation(np.radians(data['raICRS']),
                                  np.radians(data['decICRS']),
                                  ra_icrs, dec_icrs)

        dist = arcsecFromRadians(dist)
        self.assertLess(dist.max(), 0.001)

        # verify that the distance between raPhoSim, decPhoSim and
        # raICRS, decICRS is greater than the distance between
        # the original raICRS, decICRS and the newly-calculated
        # raICRS, decICRS
        dist_bad = _angularSeparation(ra_pho_rad, dec_pho_rad,
                                      np.radians(data['raICRS']),
                                      np.radians(data['decICRS']))

        dist_bad = arcsecFromRadians(dist_bad)
        self.assertGreater(dist_bad.min(), dist.max())

        del db

    def test_stellar_astrometry_degrees(self):
        """
        Test that we can go from raPhoSim, decPhoSim to ICRS coordinates
        in the case of stars (in degrees)
        """
        db = testStarsDBObj(driver='sqlite', database=self.db_name)
        cat = StarTestCatalog(db, obs_metadata=self.obs)
        with lsst.utils.tests.getTempFilePath('.txt') as cat_name:
            cat.write_catalog(cat_name)
            dtype = np.dtype([('raICRS', float), ('decICRS', float),
                             ('raPhoSim', float), ('decPhoSim', float),
                             ('raJ2000', float), ('decJ2000', float),
                             ('pmRA', float), ('pmDec', float),
                             ('parallax', float), ('vRad', float)])
            data = np.genfromtxt(cat_name, dtype=dtype)
        self.assertGreater(len(data), 100)

        # verify that, when transforming back to ICRS, we are within
        # 10^-3 arcsec
        ra_icrs, dec_icrs = PhoSimAstrometryBase.icrsFromPhoSim(data['raPhoSim'],
                                                                data['decPhoSim'],
                                                                self.obs)
        dist = angularSeparation(data['raICRS'], data['decICRS'],
                                 ra_icrs, dec_icrs)

        dist = 3600.0*dist
        self.assertLess(dist.max(), 0.001)

        # verify that the distance between raPhoSim, decPhoSim and
        # raICRS, decICRS is greater than the distance between
        # the original raICRS, decICRS and the newly-calculated
        # raICRS, decICRS
        dist_bad = angularSeparation(data['raPhoSim'], data['decPhoSim'],
                                     data['raICRS'], data['decICRS'])

        dist_bad = 3600.0*dist_bad
        self.assertGreater(dist_bad.min(), dist.max())

        del db

    def test_stellar_observed_radians(self):
        """
        Test ability to go all the way to observed RA, Dec
        from PhoSim (this is necessary for the ImSim software
        that DESC is working on)
        """
        db = testStarsDBObj(driver='sqlite', database=self.db_name)
        cat = StarTestCatalog(db, obs_metadata=self.obs)
        with lsst.utils.tests.getTempFilePath('.txt') as cat_name:
            cat.write_catalog(cat_name)
            dtype = np.dtype([('raICRS', float), ('decICRS', float),
                             ('raPhoSim', float), ('decPhoSim', float),
                             ('raJ2000', float), ('decJ2000', float),
                             ('pmRA', float), ('pmDec', float),
                             ('parallax', float), ('vRad', float)])
            data = np.genfromtxt(cat_name, dtype=dtype)
        self.assertGreater(len(data), 100)

        (ra_obs,
         dec_obs) = _observedFromICRS(data['raJ2000'],
                                      data['decJ2000'],
                                      obs_metadata=self.obs,
                                      pm_ra=data['pmRA'],
                                      pm_dec=data['pmDec'],
                                      parallax=data['parallax'],
                                      v_rad=data['vRad'],
                                      includeRefraction=True,
                                      epoch=2000.0)

        (ra_appGeo,
         dec_appGeo) = PhoSimAstrometryBase._appGeoFromPhoSim(np.radians(data['raPhoSim']),
                                                              np.radians(data['decPhoSim']),
                                                              self.obs)

        (ra_obs_2,
         dec_obs_2) = _observedFromAppGeo(ra_appGeo, dec_appGeo,
                                          obs_metadata=self.obs,
                                          includeRefraction=True)

        np.testing.assert_array_almost_equal(ra_obs, ra_obs_2, decimal=10)
        np.testing.assert_array_almost_equal(dec_obs, dec_obs_2, decimal=10)

    def test_stellar_observed_degrees(self):
        """
        Test ability to go all the way to observed RA, Dec
        from PhoSim (this is necessary for the ImSim software
        that DESC is working on)
        """
        db = testStarsDBObj(driver='sqlite', database=self.db_name)
        cat = StarTestCatalog(db, obs_metadata=self.obs)
        with lsst.utils.tests.getTempFilePath('.txt') as cat_name:
            cat.write_catalog(cat_name)
            dtype = np.dtype([('raICRS', float), ('decICRS', float),
                             ('raPhoSim', float), ('decPhoSim', float),
                             ('raJ2000', float), ('decJ2000', float),
                             ('pmRA', float), ('pmDec', float),
                             ('parallax', float), ('vRad', float)])
            data = np.genfromtxt(cat_name, dtype=dtype)
        self.assertGreater(len(data), 100)

        (ra_obs,
         dec_obs) = observedFromICRS(np.degrees(data['raJ2000']),
                                     np.degrees(data['decJ2000']),
                                     obs_metadata=self.obs,
                                     pm_ra=arcsecFromRadians(data['pmRA']),
                                     pm_dec=arcsecFromRadians(data['pmDec']),
                                     parallax=arcsecFromRadians(data['parallax']),
                                     v_rad=data['vRad'],
                                     includeRefraction=True,
                                     epoch=2000.0)

        (ra_appGeo,
         dec_appGeo) = PhoSimAstrometryBase.appGeoFromPhoSim(data['raPhoSim'],
                                                             data['decPhoSim'],
                                                             self.obs)

        (ra_obs_2,
         dec_obs_2) = observedFromAppGeo(ra_appGeo, dec_appGeo,
                                         obs_metadata=self.obs,
                                         includeRefraction=True)

        np.testing.assert_array_almost_equal(ra_obs, ra_obs_2, decimal=8)
        np.testing.assert_array_almost_equal(dec_obs, dec_obs_2, decimal=8)


    def test_galaxy_astrometry_radians(self):
        """
        Test that we can go from raPhoSim, decPhoSim to ICRS coordinates
        in the case of galaxies (in radians)
        """
        db = testGalaxyDiskDBObj(driver='sqlite', database=self.db_name)
        cat = GalaxyTestCatalog(db, obs_metadata=self.obs)
        with lsst.utils.tests.getTempFilePath('.txt') as cat_name:
            cat.write_catalog(cat_name)
            dtype = np.dtype([('raICRS', float), ('decICRS', float),
                             ('raPhoSim', float), ('decPhoSim', float)])
            data = np.genfromtxt(cat_name, dtype=dtype)
        self.assertGreater(len(data), 100)
        ra_pho_rad = np.radians(data['raPhoSim'])
        dec_pho_rad = np.radians(data['decPhoSim'])

        # verify that, when transforming back to ICRS, we are within
        # 10^-3 arcsec
        ra_icrs, dec_icrs = PhoSimAstrometryBase._icrsFromPhoSim(ra_pho_rad,
                                                                 dec_pho_rad,
                                                                 self.obs)
        dist = _angularSeparation(np.radians(data['raICRS']),
                                  np.radians(data['decICRS']),
                                  ra_icrs, dec_icrs)

        dist = arcsecFromRadians(dist)
        self.assertLess(dist.max(), 0.001)

        # verify that the distance between raPhoSim, decPhoSim and
        # raICRS, decICRS is greater than the distance between
        # the original raICRS, decICRS and the newly-calculated
        # raICRS, decICRS
        dist_bad = _angularSeparation(ra_pho_rad, dec_pho_rad,
                                      np.radians(data['raICRS']),
                                      np.radians(data['decICRS']))

        dist_bad = arcsecFromRadians(dist_bad)
        self.assertGreater(dist_bad.min(), dist.max())

        del db

    def test_galaxy_astrometry_degrees(self):
        """
        Test that we can go from raPhoSim, decPhoSim to ICRS coordinates
        in the case of galaxies (in degrees)
        """
        db = testGalaxyDiskDBObj(driver='sqlite', database=self.db_name)
        cat = GalaxyTestCatalog(db, obs_metadata=self.obs)
        with lsst.utils.tests.getTempFilePath('.txt') as cat_name:
            cat.write_catalog(cat_name)
            dtype = np.dtype([('raICRS', float), ('decICRS', float),
                             ('raPhoSim', float), ('decPhoSim', float)])
            data = np.genfromtxt(cat_name, dtype=dtype)
        self.assertGreater(len(data), 100)

        # verify that, when transforming back to ICRS, we are within
        # 10^-3 arcsec
        ra_icrs, dec_icrs = PhoSimAstrometryBase.icrsFromPhoSim(data['raPhoSim'],
                                                                data['decPhoSim'],
                                                                self.obs)
        dist = angularSeparation(data['raICRS'], data['decICRS'],
                                 ra_icrs, dec_icrs)

        dist = 3600.0*dist
        self.assertLess(dist.max(), 0.001)

        # verify that the distance between raPhoSim, decPhoSim and
        # raICRS, decICRS is greater than the distance between
        # the original raICRS, decICRS and the newly-calculated
        # raICRS, decICRS
        dist_bad = angularSeparation(data['raPhoSim'], data['decPhoSim'],
                                     data['raICRS'], data['decICRS'])

        dist_bad = 3600.0*dist_bad
        self.assertGreater(dist_bad.min(), dist.max())

        del db

    def test_naive_focal_plane_position(self):
        """
        Test deprecession of PhoSim coordinates by comparing
        the focal plane position predicted by CatSim from ICRS
        with the focal plane position predicted by CatSim from deprecessed
        coordinates.
        """

        phosim_mixin = PhoSimAstrometryBase()

        mjd = 59587.2

        # create site with no atmosphere so that we can avoid
        # refraction
        site = Site(name="LSST", pressure=0.0, humidity=0.0)

        obs = ObservationMetaData(mjd=mjd, site=site)
        ra, dec = raDecFromAltAz(31.0, 112.0, obs)

        d_sun = distanceToSun(ra, dec, obs.mjd)
        self.assertGreater(d_sun, 45.0)

        obs = ObservationMetaData(pointingRA=ra, pointingDec=dec,
                                  rotSkyPos=27.3, mjd=mjd,
                                  site=site)
        ra_icrs = np.arange(obs.pointingRA-2.0, obs.pointingRA+2.0, 0.05)
        dec_icrs = np.arange(obs.pointingDec-2.0, obs.pointingDec+2.0, 0.05)

        coord_grid = np.meshgrid(ra_icrs, dec_icrs)
        ra_icrs = coord_grid[0].flatten()
        dec_icrs = coord_grid[1].flatten()

        (xpup_icrs,
         ypup_icrs) = pupilCoordsFromRaDec(ra_icrs, dec_icrs,
                                           obs_metadata=obs,
                                           epoch=2000.0,
                                           includeRefraction=False)

        (x_focal_icrs,
         y_focal_icrs) = focalPlaneCoordsFromPupilCoords(xpup_icrs,
                                                         ypup_icrs,
                                                         camera=lsst_camera())

        ra_obs, dec_obs = observedFromICRS(ra_icrs, dec_icrs, obs_metadata=obs,
                                           epoch=2000.0,
                                           includeRefraction=False)

        ra_obs_rad = np.radians(ra_obs)
        dec_obs_rad = np.radians(dec_obs)

        (ra_deprecessed_rad,
         dec_deprecessed_rad) = phosim_mixin._dePrecess(ra_obs_rad,
                                                        dec_obs_rad, obs)

        (xpup_deprecessed,
         ypup_deprecessed) = _naivePupilCoordsFromObserved(ra_deprecessed_rad,
                                                           dec_deprecessed_rad,
                                                           obs._pointingRA,
                                                           obs._pointingDec,
                                                           obs._rotSkyPos)

        (x_focal_deprecessed,
         y_focal_deprecessed) = focalPlaneCoordsFromPupilCoords(xpup_deprecessed,
                                                                ypup_deprecessed,
                                                                camera=lsst_camera())

        dd = np.sqrt((x_focal_icrs-x_focal_deprecessed)**2
                     +(y_focal_icrs-y_focal_deprecessed)**2)

        self.assertLess(dd.max(), 1.0e-8)

    def test_against_catalog(self):
        """
        Compare deprecession results to a catalog that was validated
        with PhoSim.
        """
        obs = ObservationMetaData(pointingRA=53.00913847303155535,
                                  pointingDec=-27.43894880881512321,
                                  rotSkyPos=256.75075318193080420,
                                  mjd=59580.13955500000156462,
                                  site=Site(name="LSST", pressure=0.0,
                                            humidity=0.0))

        dtype = np.dtype([('id', int), ('ra', float), ('dec', float),
                          ('ra_deprecessed', float), ('dec_deprecessed', float),
                          ('x_dm', float), ('y_dm', float),
                          ('x_focal', float), ('y_focal', float),
                          ('x_cam', float), ('y_cam', float)])

        data = np.genfromtxt(os.path.join(getPackageDir('sims_catUtils'),
                                          'tests', 'testData',
                                          'pixel_prediction_catalog.txt'),
                             dtype=dtype)

        ra_obs, dec_obs = observedFromICRS(data['ra'], data['dec'],
                                           obs_metadata=obs,
                                           includeRefraction=False,
                                           epoch=2000.0)

        phosim_mixin = PhoSimAstrometryBase()
        ra_dep, dec_dep = phosim_mixin._dePrecess(np.radians(ra_obs),
                                                  np.radians(dec_obs),
                                                  obs)

        np.testing.assert_array_almost_equal(data['ra_deprecessed'],
                                             np.degrees(ra_dep), decimal=10)

        np.testing.assert_array_almost_equal(data['dec_deprecessed'],
                                             np.degrees(dec_dep), decimal=10)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
