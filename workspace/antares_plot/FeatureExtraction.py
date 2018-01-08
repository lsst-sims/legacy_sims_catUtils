import numpy as np
from lsst.sims.coordUtils import Sed
import gatspy

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

def skewness_and_kurtosis(flux):
    mean_flux = np.mean(flux)
    stdev_flux = np.std(flux)
    skew = np.sum((flux-mean_flux)**3)/(len(flux)*stdev_flux**3)
    kurt = np.sum((flux-mean_flux)**4)/(len(flux)*stdev_flux**4)
    return skew, kurt, stdev_flux/mean_flux

def median_absolute_deviation(flux):
    return np.median(np.abs(flux-np.median(flux)))

def stetson_k(flux, sigma_flux):

    wgt = 1.0/sigma_flux
    wgt_sum = wgt.sum()
    mean = (wgt*flux).sum()/wgt_sum
    delta = (flux-mean)/sigma_flux
    n_flux = len(flux)
    delta *= np.sqrt(n_flux/(n_flux-1.0))

    k = (np.abs(delta).sum()/n_flux)/np.sqrt((delta**2).sum()/n_flux)
    return k

def von_neumann_ratio(flux):
    """
    As per von Neumann 1941
    'Distribution of theRatio of the Mean Square Successive
    Difference to the Variance'
    Annals of Mathematical Statistics 12, 367
    """
    mean_flux = np.mean(flux)
    n_flux = len(flux)
    var_flux = np.sum((flux-mean_flux)**2)/n_flux
    delsq = np.sum((flux[1:]-flux[:-1])**2)/(n_flux-1.0)
    return delsq/var_flux

def periodic_features(time, flux, sigma_flux):
    """
    Returns
    period_uncertainty
    period
    period_SNR
    False Alarm Probability
    """

    p_min = 0.1  # minimum period of 0.1 days
    p_max = 3.0*366.0
    d_p = 0.01
    period_arr = np.arange(p_min, p_max, d_p)

    ls = gatspy.periodic.LombScargleFast().fit(time,flux,sigma_flux)
    ls_p = ls.periodogram(periods=period_arr)

    best_dex = np.argmax(ls_p)
    best_period = period_arr[best_dex]
    best_power = ls_p[best_dex]

    mean_power = np.mean(ls_p)
    stdev_power = np.std(ls_p)
    d_from_mean = np.abs(ls_p-mean_power/stdev_power)
    dex_below_mean_p_1 = np.where(d_from_mean<1.0)[0]
    dex_before = dex_below_mean_p_1[np.where(dex_below_mean_p_1<best_dex)]
    d0 = dex_before.max()
    dex_after = dex_below_mean_p_1[np.where(dex_below_mean_p_1>best_dex)]
    d1 = dex_after.min()
    period_uncertainty = 0.5*(period_arr[d1]-period_arr[d0])

    median_power = np.median(ls_p)
    snr = (best_power-median_power)/stdev_power

    omega_arg = 2.0*np.pi*time/best_period
    cos_fn = np.cos(omega_arg)
    sin_fn = np.sin(omega_arg)
    sig2 = sigma_flux**2
    bb = np.array([(flux*cos_fn/sig2).sum(),
                   (flux*sin_fn/sig2).sum(),
                   (flux/sig2).sum()])

    cc = (cos_fn**2/sig2).sum()
    cs = (cos_fn*sin_fn/sig2).sum()
    c = (cos_fn/sig2).sum()
    s = (sin_fn/sig2).sum()
    mm = np.array([[cc, cs, c],
                   [cs, ss, s],
                   [c, s, (1.0/sig2).sum()]])

    coeffs = np.linalg.solve(mm, bb)
    fit_flux = coeffs[0]*cos_fn + coeffs[1]*sin_fn + coeffs[2]
    fit_chisq = ((fit_flux-flux)**2/sig2).sum()

    wgt = 1.0/sigma2
    wgt_sum = wgt.sum()
    wgt_flux_sum = (flux*wgt).sum()
    wgt_flux_mean = wgt_flux_sum/wgt
    mean_chisq = ((wgt_flux_mean-flux)**2/sig2).sum()

    zz = 0.5*(mean_chisq-fit_chisq)   # eqn 1 of Baluev (2008) MNRAS 385, 1279

    # see paragraph under equation 10 of Baluev (2008)
    DD = np.mean(time**2) - np.mean(time)**2
    Teff = np.sqrt(4.0*np.pi*DD)
    fmax = 1.0/period_arr.min()
    W = fmax*Teff

    if W > 2.0:
        ln_pmax = -W*np.exp(-1.0*zz)*np.sqrt(zz)
        pmax = np.exp(ln_pmax)
        fap = 1.0-pmax
    else:
        fap = W*np.exp(-1.0*zz)*np.sqrt(zz)

    return period_uncertainty, best_period, snr, fap
