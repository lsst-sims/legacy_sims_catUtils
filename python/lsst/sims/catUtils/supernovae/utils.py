import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
__all__ = ['settwopanel']
def settwopanel(height_ratios=[1.0, 0.3],
                width_ratios=[1., 0.],
                padding=None,
                setdifflimits=[0.9, 1.1],
                setoffset=None,
                setgrid=[True, True],
                figsize=None):
    """
    returns a figure and axes for a main panel and a lower panel for
    showing differential information of overplotted quantities in
    the top panel.

    Parameters
    ----------
    height_ratios: list of floats, optional defaults to [1.0, 0.3]
        height ratio between the upper and lower panel
    width_ratios: list of floats, optional defaults to  [1.0, 0.0]
        width ratio between the left and right  panel
    figsize: figure size
    setgrid: List of bools, optional, defaults to [True, True]
                    whether to set grid on the two panels
    Returns
    -------
    figure object , ax0 (axes for top panel) , and ax1 (axes for lower panel)

    Examples
    --------
    >>> myfig,  myax0 , myax1 = settwopanel ( )
    >>> myax0.plot( x,  y)
    >>> myax1.plot(x, x)
    >>> myfig.tight_layout()
    """
    majorformatter = ticker.ScalarFormatter(useOffset=False)

    if figsize == None:
        fig = plt.figure()
    else:
        fig = plt.figure(figsize=figsize)

    gs = gridspec.GridSpec(
        2, 1, width_ratios=width_ratios, height_ratios=height_ratios)

    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])

    if setdifflimits != None:
        ax1.set_ylim(setdifflimits)

    ax0.set_xticklabels("", visible=False)
    ax1.yaxis.set_major_formatter(majorformatter)

    if setgrid[0]:

        ax0.grid(True)

    if setgrid[1]:
        ax1.grid(True)

    hpad = 0.0
    return fig, ax0, ax1
