import numpy as np
from lsst.sims.coordUtils import Sed

def entropy(flux, sigma_flux):
    """
    Equations 10 and 11 of Mislis et al. 2016 (MNRAS 455, 626)
    """
    if not hasattr(entropy, '_exptsq'):
        entropy._t = np.arange(0.0, 10.0, 0.01)
        exptsq = np.exp(-(entropy._t**2))
        exptsq_int = 0.5*(exptsq[:-1]+exptsq[1:])*(entropy._t[1:]-entropy._t[:-1])
        exptsq_int = np.append(np.array([0.0]), exptsq_int)
        entropy._exptsq = np.cumsum(exptsq_int)
        assert len(entropy._exptsq) == len(entropy._t)
        assert np.abs(entropy._exptsq[-1]-np.sqrt(np.pi)/2.0)<0.0001
        assert entropy._exptsq[0] == 0

    mean_flux = np.mean(flux)
    rms_flux = np.sqrt(np.sum((flux-mean_flux)**2)/(len(flux)-1.0))

    ee = 0.0
    for i_pt in range(len(flux)):
        min_flux = max(flux[i_pt] - sigma_flux[i_pt], 0.0)
        max_flux = flux[i_pt] + sigma_flux[i_pt]
        dflux = 0.01*(max_flux-min_flux)
        local_flux = np.arange(min_flux, max_flux, dflux)

        t_arg = (local_flux-mean_flux)/(rms_flux*np.sqrt(2.0))
        p_x = 0.5*(1.0 +
                   2.0*np.interp(t_arg, entropy._t, entropy._exptsq)/np.sqrt(np.pi))

        p_x_int = (0.5*(p_x[1:]+p_x[:-1])*(t_arg[1:]-t_arg[:-1])).sum()
        ee -= p_x_int
    return(ee)

def hlratio(flux):
    """
    Ratio of mean amplitude of points above the mean to mean amplitude of points
    below the mean.
    (Not sure if I should be taking those two sub-means; the description of the
    quantity in the ANTARES paper is vague)
    """
    mean_flux = np.mean(flux)
    above = np.where(flux>mean_flux)
    below = np.where(flux<mean_flux)
    mean_above = np.mean(flux[above]-mean_flux)
    mean_below = np.mean(mean_flux-flux[below])
    return mean_above/mean_below

def quartile_range(flux):
    """
    difference in magnitudes between the 25th and 75th percentile
    """
    dummy_sed = Sed()
    mag = np.sort(dummy_sed.magFromFlux(flux))
    return mag[len(mag)//4]-mag[3*len(mag)//4]
