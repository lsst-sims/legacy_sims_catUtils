import numpy as np
from lsst.sims.utils import angularSeparation
from lsst.sims.coordUtils import chipNameFromRaDecLSST

def apply_focal_plane(ra, dec, photometry_mask_1d,
                      obs_md_list, filter_obs, proper_chip):

    fov_radius = 1.75
    n_obj = len(ra)
    n_t = len(filter_obs)

    if not proper_chip:
        seed = int(np.round(ra[0]*1000.0+
                   dec[0]*1000.0+
                   obs_md_list[0].mjd.TAI))
        rng = np.random.RandomState(seed)

    chip_mask = np.zeros((n_obj, n_t), dtype=bool)
    for i_t, (obs, i_bp) in enumerate(zip(obs_md_list, filter_obs)):
        if proper_chip:
            chip_name= chipNameFromRaDecLSST(ra[photometry_mask_1d],
                                             dec[photometry_mask_1d],
                                             obs_metadata=obs,
                                             band='ugrizy'[i_bp])

            valid_chip = (np.char.find(chip_name.astype(str), 'None') == -1)
        else:
            ra_valid = ra[photometry_mask_1d]
            dec_valid = dec[photometry_mask_1d]
            dd = angularSeparation(ra_valid,
                                   dec_valid,
                                   obs.pointingRA, obs.pointingDec)

            #focal_plane_roll = rng.random_sample(len(ra_valid))
            #valid_chip = (dd<=fov_radius)&(focal_plane_roll<0.9)
            valid_chip = (dd<=fov_radius)

        chip_mask[photometry_mask_1d, i_t] = valid_chip

    return chip_mask
