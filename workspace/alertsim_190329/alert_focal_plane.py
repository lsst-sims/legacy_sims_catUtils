import numpy as np
from lsst.sims.utils import angularSeparation
from lsst.sims.coordUtils import chipNameFromRaDecLSST

def apply_focal_plane(chunk, photometry_mask_1d,
                      obs_md_list, filter_obs, proper_chip):

    fov_radius = 1.75
    n_obj = len(chunk)
    n_t = len(filter_obs)

    chip_mask = np.zeros((n_obj, n_t), dtype=bool)
    for i_t, (obs, i_bp) in enumerate(zip(obs_md_list, filter_obs)):
        if proper_chip:
            chip_name= chipNameFromRaDecLSST(chunk['ra'][photometry_mask_1d],
                                             chunk['dec'][photometry_mask_1d],
                                             obs_metadata=obs,
                                             band='ugrizy'[i_bp])

            valid_chip = (np.char.find(chip_name.astype(str), 'None') == -1)
        else:
            dd = angularSeparation(chunk['ra'][photometry_mask_1d],
                                   chunk['dec'][photometry_mask_1d],
                                   obs.pointingRA, obs.pointingDec)

            valid_chip = (dd<=fov_radius)

        chip_mask[photometry_mask_1d, i_t] = valid_chip

    return chip_mask
