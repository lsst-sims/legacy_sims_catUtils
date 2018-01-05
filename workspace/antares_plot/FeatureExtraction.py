import numpy as np

__all__ = ["autocorrelation_integral"]

def autocorrelation_integral(obshistid, flux):
    """
    From equation (9) of Mislis et al. 2016 (MNRAS 455, 626).
    They assume a regularly sampled light curve.  I am going
    to naively adopt their expression and hope that it does the
    right thing for irregularly sampled light curves.  That would
    require testing if we actually wanted to use this method for
    research results.

    Parameters
    ----------
    obshistid -- np.array of obshistid values (time as an integer)

    flux -- np.array of flux values
    """

    bad_val = flux.min()-9999.0
    mean_flux = np.mean(flux)
    rms_sq = np.sum((flux-mean_flux)**2)/(len(flux)-1)
    min_time = obshistid.min()
    max_time = obshistid.max()
    regular_time = np.arange(0,max_time-min_time, dtype=int)
    regular_flux =  bad_val*np.ones(len(regular_time))
    regular_flux[obshistid-min_time] = flux

    a_i = 0.0
    n_t = len(regular_time)
    for tau in range(1, n_t):
        xi = regular_flux[:n_t-tau]
        xi_tau = regular_flux[tau:n_t]
        assert len(xi) == len(xi_tau)
        valid = np.where(np.logical_and(xi>bad_val+10.0, xi_tau>bad_val+10.0))
        if len(valid[0]) == 0:
            continue
        xi = xi[valid] - mean_flux
        xi_tau = xi_tau[valid] - mean_flux
        numerator = (xi*xi_tau).sum()
        a_i += numerator/(len(xi)*rms_sq)

    return a_i
