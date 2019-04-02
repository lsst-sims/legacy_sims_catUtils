import os
import lsst.sims.photUtils.SignalToNoise as SNR
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed

(lsst_bp,
 lsst_hw_bp) = BandpassDict.loadBandpassesFromFiles()

sky_sed_name = os.path.join(os.environ['THROUGHPUTS_DIR'], 'baseline',
                            'darksky.dat')

assert os.path.isfile(sky_sed_name)

sky_sed = Sed()
sky_sed.readSED_flambda(sky_sed_name)

# taken from table 1 of the overview paper;
# divide by 10 to make coadds from first year
first_yr_visits = {}
first_yr_visits['u'] = 6
first_yr_visits['g'] = 8
first_yr_visits['r'] = 18
first_yr_visits['i'] = 18
first_yr_visits['z'] = 16
first_yr_visits['y'] = 16

seeing = {}
seeing['u'] = 0.92
seeing['g'] = 0.87
seeing['r'] = 0.83
seeing['i'] = 0.80
seeing['z'] = 0.78
seeing['y'] = 0.76

coadd_m5 = {}
for bp in 'ugrizy':
    phot_params = PhotometricParameters(nexp=1,
                                        exptime=first_yr_visits[bp]*30.0)

    m5_val = SNR.calcM5(sky_sed, lsst_bp[bp], lsst_hw_bp[bp],
                        phot_params, seeing[bp])

    coadd_m5[bp] = m5_val
    print(bp,m5_val)


sky_brightness = lsst_bp.magListForSed(sky_sed)
print(sky_brightness)
