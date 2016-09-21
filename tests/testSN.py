"""
SNObject_tests:
A Class containing tests to check crictical functionality for SNObject.py

The following functionality is tested:

    - SED (flambda) for unextincted SEDs in SNCosmo and SNObject
    - SED (flambda) for MW extincted SEDs in SNCosmo and SNObject (independent
        implementations of extinction using OD94 model.)
    - Band Flux for extincted SED in r Band
    - Band Mag for extincted SED in r Band

SNIaCatalog_tests:
A Class containing tests to check crictical functionality for SNIaCatalog
"""
import os
import sqlite3
import numpy as np
import unittest

# External packages used
# import pandas as pd
from pandas.util.testing import assert_frame_equal
import sncosmo
import astropy


# Lsst Sims Dependencies
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters
from lsst.sims.photUtils import BandpassDict
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import spatiallySample_obsmetadata as sample_obsmetadata
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catalogs.db import CatalogDBObject, fileDBObject

# Routines Being Tested
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.mixins import SNIaCatalog
from lsst.sims.catUtils.utils import SNIaLightCurveGenerator

# 2016 July 28
# For some reason, the Jenkins slaves used for continuous integration
# cannot properly load the astropy config directories used by sncosmo.
# To prevent this from crashing every build, we will test whether
# the directories can be accessed and, if they cannot, use unittest.skipIf()
# to skip all of the unit tests in this file.
from astropy.config import get_config_dir

_skip_sn_tests = False
try:
    get_config_dir()
except:
    _skip_sn_tests = True


class SNObject_tests(unittest.TestCase):

    def setUp(self):
        """
        Setup tests
        SN_blank: A SNObject with no MW extinction
        """

        mydir = get_config_dir()
        print '==============================='
        print '==============================='
        print (mydir)
        print '==============================='
        print '==============================='
        # A range of wavelengths in Ang
        self.wave = np.arange(3000., 12000., 50.)
        # Equivalent wavelenths in nm
        self.wavenm = self.wave / 10.
        # Time to be used as Peak
        self.mjdobs = 571190

        # Check that we can set up a SED
        # with no extinction
        self.SN_blank = SNObject()
        self.SN_blank.setCoords(ra=30., dec=-60.)
        self.SN_blank.set(z=0.96, t0=571181, x1=2.66, c=0.353, x0=1.796e-6)
        self.SN_blank.set_MWebv(0.)

        self.SN_extincted = SNObject(ra=30., dec=-60.)
        self.SN_extincted.set(z=0.96, t0=571181, x1=2.66, c=0.353,
                              x0=1.796112e-06)

        self.SNCosmoModel = self.SN_extincted.equivalentSNCosmoModel()
        self.rectify_photParams = PhotometricParameters()
        self.lsstBandPass = BandpassDict.loadTotalBandpassesFromFiles()
        self.SNCosmoBP = sncosmo.Bandpass(wave=self.lsstBandPass['r'].wavelen,
                                          trans=self.lsstBandPass['r'].sb,
                                          wave_unit=astropy.units.Unit('nm'),
                                          name='lsst_r')

    def tearDown(self):
        del self.SNCosmoBP
        del self.SN_blank
        del self.SN_extincted

    def test_SNstatenotEmpty(self):
        """
        Check that the state of SNObject, stored in self.SNstate has valid
        entries for all keys and does not contain keys with None type Values.
        """
        myDict = self.SN_extincted.SNstate
        for key in myDict.keys():
            assert myDict[key] is not None

    def test_attributeDefaults(self):
        """
        Check the defaults and the setter properties for rectifySED and
        modelOutSideRange
        """
        snobj = SNObject(ra=30., dec=-60., source='salt2')
        self.assertEqual(snobj.rectifySED, True)
        self.assertEqual(snobj.modelOutSideTemporalRange, 'zero')

        snobj.rectifySED = False
        self.assertFalse(snobj.rectifySED, False)
        self.assertEqual(snobj.modelOutSideTemporalRange, 'zero')

    def test_raisingerror_forunimplementedmodelOutSideRange(self):
        """
        check that correct error is raised if the user tries to assign an
        un-implemented model value to
        `sims.catUtils.supernovae.SNObject.modelOutSideTemporalRange`
        """
        snobj = SNObject(ra=30., dec=-60., source='salt2')
        assert snobj.modelOutSideTemporalRange == 'zero'
        with self.assertRaises(ValueError) as context:
            snobj.modelOutSideTemporalRange = 'False'
        self.assertEqual('Model not implemented, defaulting to zero method\n',
                         context.exception.message)

    def test_rectifiedSED(self):
        """
        Check for an extreme case that the SN seds are being rectified. This is
        done by setting up an extreme case where there will be negative seds, and
        checking that this is indeed the case, and checking that they are not
        negative if rectified.
        """

        snobj = SNObject(ra=30., dec=-60., source='salt2')
        snobj.set(z=0.96, t0=self.mjdobs, x1=-3., x0=1.8e-6)
        snobj.rectifySED = False
        times = np.arange(self.mjdobs - 50., self.mjdobs + 150., 1.)
        badTimes = []
        for time in times:
            sed = snobj.SNObjectSED(time=time,
                                    bandpass=self.lsstBandPass['r'])
            if any(sed.flambda < 0.):
                badTimes.append(time)
        # Check that there are negative SEDs
        assert(len(badTimes) > 0)
        snobj.rectifySED = True
        for time in badTimes:
            sed = snobj.SNObjectSED(time=time,
                                    bandpass=self.lsstBandPass['r'])
            self.assertGreaterEqual(sed.calcADU(bandpass=self.lsstBandPass['r'],
                                                photParams=self.rectify_photParams), 0.)
            self.assertFalse(any(sed.flambda < 0.))

    def test_ComparebandFluxes2photUtils(self):
        """
        The SNObject.catsimBandFlux computation uses the sims.photUtils.sed
        band flux computation under the hood. This test makes sure that these
        definitions are in sync
        """

        snobject_r = self.SN_extincted.catsimBandFlux(
            bandpassobject=self.lsstBandPass['r'],
            time=self.mjdobs)

        # `sims.photUtils.Sed`
        sed = self.SN_extincted.SNObjectSED(time=self.mjdobs,
                                            bandpass=self.lsstBandPass['r'])
        sedflux = sed.calcFlux(bandpass=self.lsstBandPass['r'])
        np.testing.assert_allclose(snobject_r, sedflux / 3631.0)

    def test_CompareBandFluxes2SNCosmo(self):
        """
        Compare the r band flux at a particular time computed in SNObject and
        SNCosmo for MW-extincted SEDs. While the underlying sed is obtained
        from SNCosmo the integration with the bandpass is an independent
        calculation in SNCosmo  and catsim
        """

        times = self.mjdobs
        catsim_r = self.SN_extincted.catsimBandFlux(
            bandpassobject=self.lsstBandPass['r'],
            time=times)
        sncosmo_r = self.SNCosmoModel.bandflux(band=self.SNCosmoBP,
                                               time=times, zpsys='ab',
                                               zp=0.)
        np.testing.assert_allclose(sncosmo_r, catsim_r)

    def test_CompareBandMags2SNCosmo(self):
        """
        Compare the r band flux at a particular time computed in SNObject and
        SNCosmo for MW-extincted SEDs. Should work whenever the flux comparison
        above works.
        """
        times = self.mjdobs
        catsim_r = self.SN_extincted.catsimBandMag(
            bandpassobject=self.lsstBandPass['r'],
            time=times)
        sncosmo_r = self.SNCosmoModel.bandmag(band=self.SNCosmoBP,
                                              time=times, magsys='ab')
        np.testing.assert_allclose(sncosmo_r, catsim_r)

    def test_CompareExtinctedSED2SNCosmo(self):
        """
        Compare the extincted SEDS in SNCosmo and SNObject. Slightly more
        non-trivial than comparing unextincted SEDS, as the extinction in
        SNObject uses different code from SNCosmo. However, this is still
        using the same values of MWEBV, rather than reading it off a map.
        """
        SNObjectSED = self.SN_extincted.SNObjectSED(time=self.mjdobs,
                                                    wavelen=self.wavenm)

        SNCosmoSED = self.SNCosmoModel.flux(time=self.mjdobs, wave=self.wave) \
            * 10.

        np.testing.assert_allclose(SNObjectSED.flambda, SNCosmoSED,
                                   rtol=1.0e-7)

    def test_CompareUnextinctedSED2SNCosmo(self):
        """
        Compares the unextincted flux Densities in SNCosmo and SNObject. This
        is mereley a sanity check as SNObject uses SNCosmo under the hood.
        """

        SNCosmoFluxDensity = self.SN_blank.flux(wave=self.wave,
                                                time=self.mjdobs) * 10.

        unextincted_sed = self.SN_blank.SNObjectSED(time=self.mjdobs,
                                                    wavelen=self.wavenm)

        SNObjectFluxDensity = unextincted_sed.flambda
        np.testing.assert_allclose(SNCosmoFluxDensity, SNObjectFluxDensity,
                                   rtol=1.0e-7)

    def test_redshift(self):
        """
        test that the redshift method works as expected by checking that
        if we redshift a SN from its original redshift orig_z to new_z where
        new_z is smaller (larger) than orig_z:
        - 1. x0 increases (decreases)
        - 2. source peak absolute magnitude in BesselB band stays the same
        """
        from astropy.cosmology import FlatLambdaCDM
        cosmo = FlatLambdaCDM(H0=70., Om0=0.3)

        orig_z = self.SN_extincted.get('z')
        orig_x0 = self.SN_extincted.get('x0')
        peakabsMag = self.SN_extincted.source_peakabsmag('BessellB', 'AB', cosmo=cosmo)

        lowz = orig_z * 0.5
        highz = orig_z * 2.0

        # Test Case for lower redshift
        self.SN_extincted.redshift(z=lowz, cosmo=cosmo)
        low_x0 = self.SN_extincted.get('x0')
        lowPeakAbsMag = self.SN_extincted.source_peakabsmag('BessellB', 'AB', cosmo=cosmo)

        # Test 1.
        self.assertGreater(low_x0, orig_x0)
        # Test 2.
        self.assertEqual(peakabsMag, lowPeakAbsMag)

        # Test Case for higher redshift
        self.SN_extincted.redshift(z=highz, cosmo=cosmo)
        high_x0 = self.SN_extincted.get('x0')
        HiPeakAbsMag = self.SN_extincted.source_peakabsmag('BessellB', 'AB', cosmo=cosmo)

        # Test 1.
        self.assertLess(high_x0, orig_x0)
        # Test 2.
        self.assertEqual(peakabsMag, HiPeakAbsMag)

