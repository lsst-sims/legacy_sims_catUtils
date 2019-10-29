import os
import h5py
import numpy as np
import lsst.sims.photUtils.SignalToNoise as SNR
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed

out_name = 'data/dflux_SNR5_lookup.h5'

(lsst_bp,
 lsst_hw_bp) = BandpassDict.loadBandpassesFromFiles()

sky_sed_name = os.path.join(os.environ['THROUGHPUTS_DIR'], 'baseline',
                            'darksky.dat')

assert os.path.isfile(sky_sed_name)

sky_sed = Sed()
sky_sed.readSED_flambda(sky_sed_name)

# from the overview paper
#m5_single = {}
#m5_single['u'] = 23.57
#m5_single['g'] = 24.65
#m5_single['r'] = 24.21
#m5_single['i'] = 23.79
#m5_single['z'] = 23.21
#m5_single['y'] = 22.31

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


m5_single = 26.0
mag_grid = np.arange(10.0, 33.0, 0.05)
dummy_sed = Sed()

phot_params_single = PhotometricParameters(nexp=1, exptime=30.0)

with h5py.File(out_name, 'w') as out_file:
    out_file.create_dataset('mag_grid', data=mag_grid)
    for bp in 'ugrizy':
        phot_params_coadd = PhotometricParameters(nexp=1,
                                                  exptime=first_yr_visits[bp]*30.0)

        gamma_coadd = None
        gamma_single = None
        dflux_grid = []
        for mag in mag_grid:
            flux = dummy_sed.fluxFromMag(mag)

            (snr_coadd,
             gamma_coadd) = SNR.calcSNR_m5(mag, lsst_bp[bp], coadd_m5[bp],
                                           phot_params_coadd,
                                           gamma=gamma_coadd)

            noise_coadd = flux/snr_coadd

            (snr_single,
             gamma_single) = SNR.calcSNR_m5(mag, lsst_bp[bp], m5_single,
                                            phot_params_single,
                                            gamma=gamma_single)

            noise_single = flux/snr_single
            #print(noise_coadd, noise_single, noise_single/noise_coadd, np.sqrt(first_yr_visits[bp]))
            noise_tot = np.sqrt(noise_single**2+noise_coadd**2)
            dflux = 5.0*noise_tot
            dflux_grid.append(dflux)
        dflux_grid = np.array(dflux_grid)
        out_file.create_dataset('%s_dflux' % bp, data=dflux_grid)
