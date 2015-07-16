from lsst.sims.photUtils import Bandpass, Sed, PhotometricParameters

__all__ = ["calcADUwrapper"]


def calcADUwrapper(sedName=None, magNorm=None, redshift=None, internalAv=None, internalRv=None,
                   galacticAv=None, galacticRv=None, bandpass=None):

    imsimband = Bandpass()
    imsimband.imsimBandpass()
    sed = Sed()
    sed.readSED_flambda(sedName)
    fNorm = sed.calcFluxNorm(magNorm, imsimband)
    sed.multiplyFluxNorm(fNorm)
    if internalAv is not None and internalRv is not None:
        if internalAv != 0.0 and internalRv != 0.0:
            a_int, b_int = sed.setupCCMab()
            sed.addCCMDust(a_int, b_int, A_v=internalAv, R_v=internalRv)
    
    if redshift is not None and redshift != 0.0:
        sed.redshiftSED(redshift, dimming=True)
    
    a_int, b_int = sed.setupCCMab()
    sed.addCCMDust(a_int, b_int, A_v=galacticAv, R_v=galacticRv)
    
    adu = sed.calcADU(bandpass, photParams=PhotometricParameters())
    
    return adu
