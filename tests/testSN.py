#!/usr/bin/env python
"""
A Class containing tests to check crictical functionality for SNObject.py

The following functionality is tested:

    - SED (flambda) for unextincted SEDs in SNCosmo and SNObject
    - SED (flambda) for MW extincted SEDs in SNCosmo and SNObject (independent
        implementations of extinction using OD94 model.)
    - Band Flux for extincted SED in r Band
    - Band Mag for extincted SED in r Band
"""
import numpy as np

import unittest
import lsst.utils.tests as utilsTests
from lsst.sims.photUtils import Bandpass
from lsst.sims.photUtils import BandpassDict

# Routines Being Tested
from lsst.sims.catUtils.mixins import SNObject
# External package
import sncosmo
import astropy


class SNObject_tests(unittest.TestCase):

    def setUp(self):
        """
        Setup tests

        
        SN_blank: A SNObject with no MW extinction


        """

        # A range of wavelengths in Ang
        self.wave = np.arange(3000., 12000., 50.)
        # Equivalent wavelenths in nm
        self.wavenm = self.wave/10.
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
        
        self.lsstBandPass = BandpassDict.loadTotalBandpassesFromFiles()
        self.SNCosmoBP = sncosmo.Bandpass(wave=self.lsstBandPass['r'].wavelen,
                                          trans=self.lsstBandPass['r'].sb,
                                          wave_unit=astropy.units.Unit('nm'),
                                          name='lsst_r')


    def tearDown(self):
        #if os.path.exists(self.catName):
        #    os.unlink(self.catName)
        #del self.catName
        pass

    def test_bandFluxes(self):
        """
        Compare the r band flux at a particular time computed in SNObject and
        SNCosmo for MW-extincted SEDs
        """
        times = self.mjdobs
        catsim_r = self.SN_extincted.catsimBandFluxes(
                                        bandpassobject=self.lsstBandPass['r'],
                                        time=times)
        sncosmo_r = self.SNCosmoModel.bandflux(band=self.SNCosmoBP,
                                               time=times,  zpsys='ab',
                                               zp=0.)
        np.testing.assert_allclose(sncosmo_r, catsim_r)

    def test_bandMags(self):
        """
        Compare the r band flux at a particular time computed in SNObject and
        SNCosmo for MW-extincted SEDs
        """
        times = self.mjdobs
        catsim_r = self.SN_extincted.catsimBandMags(
                                        bandpassobject=self.lsstBandPass['r'],
                                        time=times)
        sncosmo_r = self.SNCosmoModel.bandmag(band=self.SNCosmoBP,
                                               time=times,  magsys='ab')
        np.testing.assert_allclose(sncosmo_r, catsim_r)
        
    def test_extinctedSED(self):
        """
        Compare the extincted SEDS in SNCosmo and SNObject. Slightly more
        non-trivial than comparing unextincted SEDS, as the extinction in
        SNObject uses different code from SNCosmo. However, this is still
        using the same values of MWEBV, rather than reading it off a map.
        """
        SNObjectSED = self.SN_extincted.SNObjectSED(time=self.mjdobs,
                                                    wavelen=self.wavenm)

        SNCosmoSED = self.SNCosmoModel.flux(time=self.mjdobs, wave=self.wave) * 10.

        np.testing.assert_allclose(SNObjectSED.flambda, SNCosmoSED, rtol=1.0e-7)

    def test_unextinctedSED(self):
        """
        Compares the unextincted flux Densities in SNCosmo and SNObject. This is mereley
        a sanity check as SNObject uses SNCosmo under the hood.
        """
        SNCosmoFluxDensity = self.SN_blank.flux(wave=self.wave,
                                                time=self.mjdobs) * 10.
        unextincted_sed = self.SN_blank.SNObjectSED(time=self.mjdobs,
                                                    wavelen=self.wavenm)
        SNObjectFluxDensity = unextincted_sed.flambda 
        np.testing.assert_allclose(SNCosmoFluxDensity, SNObjectFluxDensity,
                                   rtol=1.0e-7)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(SNObject_tests)
    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(),shouldExit)

if __name__ =='__main__':
    run(True)

