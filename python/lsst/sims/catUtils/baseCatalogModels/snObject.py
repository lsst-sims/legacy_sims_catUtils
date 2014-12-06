"""
Class describing the SN object itself. The SN object derives from
SNCosmo.Model and has additional properties such as ra, dec.
It also has additional methods :
-     calc_mags which use the magnitude calculations in LSST stack
-     extinction which use the extinction calculations in LSST stack

- Usage:  (after setups described in the readme.rst) python snObject.py

"""
import sncosmo
import numpy as np
import os
import matplotlib.pyplot as plt

from astropy.units import Unit
from astropy.coordinates import SkyCoord
from sncosmo import Model

from lsst.sims.photUtils.Photometry import PhotometryStars, Sed, Bandpass
from lsst.sims.photUtils.EBV import EBVbase

dustmaproot = os.getenv('SIMS_DUSTMAPS_DIR')
map_dir = os.path.join(dustmaproot, 'DustMaps')

wavelenstep = 0.1
plot = False


def getlsstbandpassobjs(loadsncosmo=True, loadcatsim=True):
    """
    General utility to return a list of the baseline LSST bandpasses loaded
    as catsim bandpass objects, and register them as SNCosmo bandpasses
    accessible through strings like 'LSSTu'.
    args:
        loadsncosmo: Bool, optional, defaults to True
            variable to decide whether to register the LSST bandpasses as
            SNCosmo registered bandpass objects accessible through strings
            like 'LSSTu'
        loadcatsim : Bool, optional, defaults to True
            variable to decide whether to set up catsim bandpass objects
            are return the list of u,g,r,i,z,y bandpasses
    returns:
        if loadcatsim is true, list of catsim bandpass objects corresponding
        to LSST baseline u, g, r, i, z, y filters.
        if loadcatsim is False, return is None

    Examples:

    """
    bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
    banddir = os.path.join(os.getenv('THROUGHPUTS_DIR'), 'baseline')
    lsstbands = []
    lsstbp = {}

    for band in bandPassList:
        # setup sncosmo bandpasses
        bandfname = banddir + "/total_" + band + '.dat'
        if loadsncosmo:
            # register the LSST bands to the SNCosmo registry
            # Not needed for LSST, but useful to compare independent codes
            # Usually the next two lines can be merged,
            # but there is an astropy bug currently which affects only OSX.
            numpyband = np.loadtxt(bandfname)
            sncosmoband = sncosmo.Bandpass(wave=numpyband[:, 0],
                                           trans=numpyband[:, 1],
                                           wave_unit=Unit('nm'),
                                           name='LSST' + band)
            sncosmo.registry.register(sncosmoband, force=True)
        if loadcatsim:
            # Now load LSST bandpasses for catsim
            lsstbp[band] = Bandpass()
            lsstbp[band].readThroughput(bandfname, wavelen_step=wavelenstep)
            lsstbands.append(lsstbp[band])
    fifilterfigs, filterax = plt.subplots()
    if plot:
        for band in bandPassList:
            b = sncosmo.get_bandpass('LSST' + band)
            filterax.plot(b.wave, b.trans, '-k', lw=2.0)
            filterax.set_xlabel(r'$\lambda$, ($\AA$)')
            filterax.set_ylabel(r'transmission')
        plt.show()

    if loadcatsim:
        return lsstbands
    else:
        return None


class SNObject (Model):
    """
    Extension of the SNCosmo TimeSeries Model to include more parameters and
    use methods in the catsim stack. We constrain ourselves to the use of a
    specific SALT model for the Supernova (Salt2-Extended), and set its MW
    extinction to be 0.

    new attributes:
    ---------------
    ra:
        ra of the SN in degrees
    dec:
        dec of the SN in degrees

    new methods:
    ------------
    mwebvfrommaps: Uses the LSST stack to obtain MW extinction according to
        CCM 89, the ra and dec of the supernova, and the SFD dustmaps to apply
        appropriate extinction to the SN sed. must be run after the ra, dec
        parameters are set.
        args:
        returns:
    set_mwebv(values): Set the value of attribute _mwebv to a particular
        value
    lsstbandmags: Uses the LSST stack functionality to obtain LSST band
        magnitudes using the bandpass filters.
        args:
        returns:
    Notes:
    ------
    """
    def __init__(self, ra=None, dec=None):
        Model.__init__(self, source="salt2-extended",
                       effects=[sncosmo.CCM89Dust()], effect_names=['mw'],
                       effect_frames=['obs'])
        # Current implementation of Model has a default value of mwebv = 0.
        # ie. no extinction, but this is not part of the API, so should not
        # depend on it, set explicitly in order to unextincted SED from
        # SNCosmo. We will use catsim extinction from photUtils.
        self.set(mwebv=0.)
        # If we know ra, dec in degree
        self.ra = ra
        self.dec = dec
        self.lsstmwebv = EBVbase()
        self.mwebvfrommaps()
        self._seed = None
        return

    # Comment to self: Don't see why having the seed is of any help here
    # will probably remove
    @property
    def seed(self):
        return self._seed

    def set_mwebv(self, value):
        """
        if mwebv value is known, this can be used to set the attribute
        _mwebv of the SNObject class to the value (float).
        args:
            value: float, mandatory
                value of mw extinction parameter E(B-V) in mags to be used in
                applying extinction to the SNObject spectrum
        returns:
            None
        Notes:
            For a large set of SN, one may use fast `np.ndarray` valued
            functions to obtain an array of such values, and then set the
            values from such an array.
        """
        self._mwebv = value
        return

    def mwebvfrommaps(self):
        """
        set the attribute _mwebv of the class from the ra and dec
        of the SN. If the ra or dec attribute of the class is None,
        set this attribute to None.
        args: None
        returns: None
        Notes:
            This function must be run after the class has attributes ra and
            dec specified. In case it is run before this, the mwebv value will
            be set to None.

        """
        ra = self.ra
        dec = self.dec
        if ra is None or dec is None:
            self._mwebv = None
            return
        skycoord = np.array([[ra], [dec]]) * np.pi / 180.
        # skycoords = SkyCoord(ra, dec, unit='deg')
        # t_mwebv = sncosmo.get_ebv_from_map(skycoords,
        #                                    mapdir=map_dir,
        #                                    interpolate=False)
        self._mwebv = self.lsstmwebv.calculateEbv(equatorialCoordinates=
                                                  skycoord)[0]
        # print "compare vals :", t_mwebv, self._mwebv
        return

    def lsstbandmags(self, lsstbands, time):
        """
        return a numpy array of magnitudes of the SN spectrum in the ab
        magnitude system.
        args:
            lsstbands: mandatory, list of bandpass objects
                a list of LSST catsim bandpass objects
            time: mandatory, float
                MJD at which this is evaluated
        returns:
            `np.ndarray` of mag values for each band in lsstbandpass.
            Unphysical values of the flux density are reported as np.nan
        """
        filterwav = lsstbands[0].wavelen
        SEDfromSNcosmo = Sed(wavelen=filterwav,
                             flambda=self.flux(time=time,
                                               wave=filterwav*10.)*10.)
        # Check if I need this sync
        SEDfromSNcosmo.synchronizeSED(wavelen_min=filterwav[0],
                                      wavelen_max=filterwav[-2],
                                      wavelen_step=wavelenstep)
        # Apply LSST extinction
        ax, bx = SEDfromSNcosmo.setupCCMab()
        SEDfromSNcosmo.addCCMDust(a_x=ax, b_x=bx, ebv=self._mwebv)
        phiarray, dlambda = SEDfromSNcosmo.setupPhiArray(lsstbands)
        SEDfromSNcosmo.synchronizeSED(wavelen_min=filterwav[0],
                                      wavelen_max=filterwav[-2],
                                      wavelen_step=wavelenstep)

        return SEDfromSNcosmo.manyMagCalc(phiarray, wavelen_step=wavelenstep)


if __name__ == "__main__":
    """
    Example code for writing a light curve in multiple bands.
    """

    import numpy as np
    import matplotlib.pyplot as plt

    ra = 204.
    dec = -30.
    SN = SNObject(ra, dec)
    SN.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)
    print SN
    SNCosmoSN = SNObject(ra, dec)
    SNCosmoSN.set(x0=1.847e-6, x1=0.1, z=0.2, mwebv=SN._mwebv)
    lsstbands = getlsstbandpassobjs()
    sncosmobands = ['LSSTu', 'LSSTg', 'LSSTr', 'LSSTi', 'LSSTz', 'LSSTy']
    w = sncosmo.get_bandpass(sncosmobands[0]).wave
    l = []
    for time in np.arange(-20., 50., 1.0):
        t = time*np.ones(len(sncosmobands))
        t.tolist()
        x = SN.lsstbandmags(lsstbands, time=time)
        y = SNCosmoSN.bandmag(band=sncosmobands, time=t, magsys='ab')
        e = [time]
        e += x.tolist()
        e += y.tolist()
        l.append(e)
    header = "time(mjd) u g r i z y su sg sr si sz sy"
    np.savetxt('../out/lc.dat', np.array(l), fmt='%10.6f', header=header)
