import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.catalogs.measures.instance.fileMaps import defaultSpecMap
import matplotlib.pylab as plt
import warnings
import numpy
import copy
import os

class ValidationUtils(object):
    sedDict = {}
    imsimband = Bandpass()
    imsimband.imsimBandpass()

    @classmethod
    def calcMags(cls, sedKey, magNorm, bplist, z, matchwavelen=None, internalAv=None, 
            phiarray=None, dlambda=None):
        if sedKey is None or sedKey == 'None' or magNorm is None or numpy.isnan(magNorm):
            return None, None, None, None
        datadir = os.path.join(os.environ.get("CAT_SHARE_DATA"), "data")
        if cls.sedDict.has_key(sedKey):
            sed = copy.deepcopy(cls.sedDict[sedKey])
        else:
            sed = Sed()
            try:
                sed.readSED_flambda(os.path.join(datadir, defaultSpecMap[sedKey]))
            except KeyError, ke:
                warnings.warn("Couldn't find the spectrum for one of the components")
                return None, None, None, None
            if matchwavelen is not None:
                sed.resampleSED(wavelen_match=matchwavelen)
            cls.sedDict[sedKey] = copy.deepcopy(sed)
        if phiarray is None or dlambda is None:
            phiarray, dlambda = sed.setupPhiArray(bplist)
        matchwavelen = sed.wavelen
        try:
            fNorm = sed.calcFluxNorm(magNorm, cls.imsimband)
        except Exception, e:
            warnings.warn("Could not calculate the flux norm")
            return None, None, None, None
        sed.multiplyFluxNorm(fNorm)
        #I guess this assumes rv=3.1??
        if internalAv:
            a_int, b_int = sed.setupCCMab()
            sed.addCCMDust(a_int, b_int, A_v=internalAv)
        sed.redshiftSED(z, dimming=True)
        try:
            sed.resampleSED(wavelen_match=bplist[0].wavelen)
        except Exception, e:
            warnings.warn("Exception trying to resample the SED.  Redshift is: %f"%(z))
            return None, None, None, None
        sed.flambdaTofnu()
        return sed.manyMagCalc(phiarray, dlambda), matchwavelen, phiarray, dlambda

    @staticmethod
    def sumMags(*args):
        retmags = None
        m_o = 22.
        first = True
        for arg in args:
            if arg is None:
                continue
            if first:
                retmags = numpy.power(10, (arg - m_o)/-2.5)
                first = False
            else:
                retmags += numpy.power(10, (arg - m_o)/-2.5)
        if retmags is not None:
            return -2.5*numpy.log10(retmags) + m_o
        else:
            return None

    @classmethod
    def calcTotalMags(cls, cat, bands):
        # Set up phi, the wavelength-normalized system response for each filter,
        # for each bandpass for manyMagCalc method.
        bplist = []
        for f in bands:
            bands[f].sbTophi()
            bplist.append(bands[f])

        diskfile = cat.column_by_name('sedFilenameDisk')
        bulgefile = cat.column_by_name('sedFilenameBulge')
        agnfile = cat.column_by_name('sedFilenameAgn')

        diskmn = cat.column_by_name('magNormDisk')
        bulgemn = cat.column_by_name('magNormBulge')
        agnmn = cat.column_by_name('magNormAgn')

        bulgeAv = cat.column_by_name('internalAvBulge')
        diskAv = cat.column_by_name('internalAvDisk')

        diskData = numpy.asarray(zip(diskfile, diskmn, diskAv), dtype=[('sedfile', (str, 50)), ('norm', float), ('av', float)])
        bulgeData = numpy.asarray(zip(bulgefile, bulgemn, bulgeAv), dtype=[('sedfile', (str, 50)), ('norm', float), ('av', float)])
        agnData = numpy.asarray(zip(agnfile, agnmn), dtype=[('sedfile', (str, 50)), ('norm', float)])
        redshift = cat.column_by_name('redshift')

        retMags = dict([(k, []) for k in bands])
        retbMags = dict([(k, []) for k in bands])
        retdMags = dict([(k, []) for k in bands])
        retaMags = dict([(k, []) for k in bands])
        for ddat, bdat, adat, z in zip(diskData, bulgeData, agnData, redshift):
            dmags, matchwavelen, phiarray, dlambda = cls.calcMags(ddat['sedfile'], ddat['norm'], bplist, z, internalAv=ddat['av'])
            bmags, matchwavelen, phiarray, dlambda = cls.calcMags(bdat['sedfile'], bdat['norm'], bplist, z, internalAv=bdat['av'])
            if bmags is not None:
                amags, ased, phiarray, dlambda = cls.calcMags(adat['sedfile'], adat['norm'], bplist, z, matchwavelen=matchwavelen, phiarray=phiarray, dlambda=dlambda)
            elif dmags is not None:
                amags, ased, phiarray, dlambda = cls.calcMags(adat['sedfile'], adat['norm'], bplist, z, matchwavelen=matchwavelen, phiarray=phiarray, dlambda=dlambda)
            else:
                warnings.warn("couldn not get wavelength grid to match to")
            mags = cls.sumMags(dmags, bmags, amags)
                
            for i,k in enumerate(bands):
                if mags is not None:
                    retMags[k].append(mags[i])
                else:
                    retMags[k].append(None)
                if dmags is not None:
                    retdMags[k].append(dmags[i])
                else:
                    retdMags[k].append(None)
                if bmags is not None:
                    retbMags[k].append(bmags[i])
                else:
                    retbMags[k].append(None)
                if amags is not None:
                    retaMags[k].append(amags[i])
                else:
                    retaMags[k].append(None)
        return retMags, retbMags, retdMags, retaMags
    
    @classmethod
    def calcTotalLSSTMags(cls, cat, bandpasses=('u', 'g', 'r', 'i', 'z', 'y')):
        tpath = os.getenv('LSST_THROUGHPUTS_DEFAULT')
        bands = {}
        for k in bandpasses:
            bands[k] = Bandpass()
            bands[k].readThroughput(os.path.join(tpath, "total_%s.dat"%k))
        return cls.calcTotalMags(cat, bands)

    @classmethod
    def calcTotalSDSSMags(cls, cat, bandpasses=('u', 'g', 'r', 'i', 'z')):
        tpath = os.getenv('SDSS_THROUGHPUTS')
        bands = {"u":None, "g":None, "r":None, "i":None, "z":None}
        for k in bandpasses:
            bands[k] = Bandpass()
            bands[k].readThroughput(os.path.join(tpath, "sdss_%s.dat"%k))
        return cls.calcTotalMags(cat, bands)

    @staticmethod
    def read_Durham_counts(filename):
        fh = open(filename)
        datax = []
        datay = []
        data_errors = []
        errp = []
        errn = []
        for l in fh:
            if l.startswith("#"):
                continue
            els = l.rstrip().split()
            if len(els) == 0:
                continue
            #There are several cases following are the assumptions on len(els):
            #2: x, y
            #3: x, y, y_corr
            #4: x, y, -ve_err, +ve_err
            #5: x, y, y_corr, -ve_err, +ve_err
            if len(els) == 2:
                continue
                x = float(els[0])
                y = numpy.power(10.,float(els[1]))
                err1 = numpy.power(10.,0.)
                err2 = numpy.power(10.,0.)
            elif len(els) == 3:
                continue
                x = float(els[0])
                y = numpy.power(10.,float(els[2]))
                err1 = numpy.power(10.,0.)
                err2 = numpy.power(10.,0.)
            elif len(els) == 4:
                x = float(els[0])
                y = numpy.power(10.,float(els[1]))
                err1 = numpy.power(10.,float(els[3]))
                err2 = numpy.power(10.,float(els[2]))
            elif len(els) == 5:
                x = float(els[0])
                y = numpy.power(10.,float(els[2]))
                err1 = numpy.power(10.,float(els[4]))
                err2 = numpy.power(10.,float(els[3]))
            else:
                raise ValueError("Not one of the expected configurations")
            datax.append(x)
            datay.append(y)
            errn.append(y - y/err1)
            errp.append(y*err2-y)
        data_errors=[errn, errp]
        return datax, datay, data_errors
        
    @staticmethod
    def plotN(compx, compy, comperr, basemean, binx, filename=None):
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(111)
        dfit = numpy.polyfit(compx, numpy.log10(compy), 2)
        ax.errorbar(compx, compy, yerr=comperr, fmt='or', alpha=0.5, label='Empirical Data')
        ax.plot(binx, numpy.power(10,numpy.polyval(dfit, binx)), c='k', label='Quadratic fit')
        ax.plot(binx, basemean, 'sb', label='Binned Base Catalog')
        ax.set_yscale('log')
        ax.axvline(x=26.8, lw=3, c='k', ls='--')
        #ax.set_grid()
        ax.set_ylim(100, 200000)
        ax.set_xlim(18., 27.)
        ax.legend(loc='upper left')
        plt.ylabel("N/0.5mag/deg$^2$")
        plt.xlabel("i-band Mag")
        if filename is not None:
    	   plt.savefig(filename)
        else:
    	   plt.show()

    @staticmethod
    def plotNResidCume(compx, compy, comperr, basemean, binx, filename=None):
        dfit = numpy.polyfit(compx, numpy.log10(compy), 2)
        fig = plt.figure(figsize=(12,9))
        ax = fig.add_subplot(111)
        comp_fit = numpy.power(10, numpy.polyval(dfit, binx))
        bins = binx - (binx[1]-binx[0])/2.
        bins = numpy.append(bins,binx[-1] + binx[1]-binx[0])
        #Here we are assuming the negative error is representative
        comp_err = numpy.histogram(compx, bins=bins, weights=comperr[0]*comperr[0])[0]
        comp_cume = numpy.cumsum(comp_fit)
        base_cume = numpy.cumsum(basemean)
        comp_cume_err=numpy.sqrt(numpy.cumsum(comp_err))/comp_cume
        ax.plot(binx, base_cume/comp_cume, 'or', label='Binned Base Catalog')
        ax.errorbar(binx, comp_cume/comp_cume, yerr=comp_cume_err, c='k', alpha=1., label="Binned Errorbars")
        ax.axvline(x=26.8, lw=3, c='k', ls='--')
        ax.axhline(y=1.1, lw=3, c='k', ls='--')
        ax.axhline(y=0.9, lw=3, c='k', ls='--')
        plt.ylim((0.7,1.3))
        plt.xlim(18., 27.)
        plt.ylabel("N$_{Base}$(<m)/N$_{Durham}$(<m)")
        plt.xlabel("i-band Mag")
        ax.legend(loc='upper right')
        if filename is not None:
    	   plt.savefig(filename)
        else:
	       plt.show()

    @staticmethod
    def intNofz(z_o, z):
        return -z_o*numpy.exp(-z/z_o)*(2.*z_o**2 + 2.*z_o*z + z**2)

    @classmethod
    def coilNofz(cls, z, z_o, norm=1., zlimits=(0.,30.)):
        area = cls.intNofz(z_o, zlimits[1]) - cls.intNofz(z_o, zlimits[0])
        return norm*z**2*numpy.exp(-z/z_o)/area

    @classmethod
    def sampleNofz(cls, z_o, n, seed=None):
        dz = 0.0001
        zs = numpy.arange(0.,6.,0.0001)
        vals = cls.coilNofz(zs, z_o)
        cumevals = numpy.cumsum(vals*dz)
        if seed is not None:
            numpy.random.seed(seed)
        rands = numpy.random.random(n)
        zsample = numpy.interp(rands, cumevals, zs)
        return zsample


    @classmethod
    def plotNofz(cls, histData, zbinx, magbins, z_o_corr, z_o_err, area, compfunc=None, filename=None):
        if compfunc is None:
            compfunc = cls.coilNofz
        fig = plt.figure(figsize=(15,15))
        barwidth = zbinx[1]-zbinx[0]
        for i in xrange(len(magbins[:-1])):
            #hist = [el[i] for el in histData]
            hist = histData[i]
            hist = numpy.array(hist)
            meanz = numpy.sum(hist*zbinx)/numpy.sum(hist)
            nsamp = 100000
            meanzsamp = numpy.sum(cls.sampleNofz(z_o_corr[i],nsamp))/nsamp
            print meanz, meanzsamp
            hist /= area
            ax=fig.add_subplot(3,2,i+1)
            plt.bar(zbinx, hist, width=barwidth, label="Base Catalog")
            norm = numpy.sum(hist)*(zbinx[1]-zbinx[0])
            cpts = numpy.arange(zbinx[0], zbinx[-1], 0.05)    
            plt.plot(cpts, compfunc(cpts, z_o_corr[i]+3*z_o_err[i], norm=norm), label="z$_o$+=3$\sigma$", ls="--", c='r', lw=2)
            plt.plot(cpts, compfunc(cpts, z_o_corr[i], norm=norm), label="z$_o$ = z$_o$", c='k', lw=2)
            plt.plot(cpts, compfunc(cpts, z_o_corr[i]-3*z_o_err[i], norm=norm), label="z$_o$ -= 3$\sigma$", ls="-.", c='r', lw=2)
            plt.ylabel("N/deg/$3$x$10^4$s$^{-1}$km")
            plt.xlabel("redshift")
            plt.title("%i < i < %i: \n base catalog mean = %f model mean = %f"%(magbins[i], magbins[i+1], meanz, meanzsamp))
            if i == 0:
                plt.legend(loc="upper right")
        plt.tight_layout()
        if filename is not None:
            plt.savefig(filename)
        else:
            plt.show()

