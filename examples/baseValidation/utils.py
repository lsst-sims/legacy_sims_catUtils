import scipy
from scipy.interpolate import UnivariateSpline
import numpy
import os
import time
import sys
import pickle
import matplotlib.pyplot as plt
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.sims.photUtils.Sed as Sed
import lsst.sims.photUtils.Bandpass as Bandpass
from lsst.sims.photUtils import EBV
import copy

def read_ACS_cat(filename, limit=None):
    fh = open(filename)
    cols = None
    data = {}
    i = 0
    for l in fh:
        if l.startswith('\\'):
            continue
        if l.startswith('|'):
            if cols is None:
                cols = [el.strip() for el in l.rstrip().split('|')[1:-1]]
                for col in cols:
                    data[col] = []
            continue
        if limit & i > limit:
	    continue
        els = l.rstrip().split()
        for col, el in zip(cols,els):
            if col == 'field':
                data[col].append(el)
            elif el == 'null':
                data[col].append(None)
            else:
                data[col].append(float(el))
        i += 1
    for col in cols:
        data[col] = numpy.asarray(data[col])
    return data
        

def read_HDF_radius(filename):
    fh = open(filename)
    mag1 = []
    mag2 = []
    rad = []
    for l in fh:
        if l.startswith("#"):
	    continue
        els = l.rstrip().split()
	if len(els) < 10:
	    continue
        mag1.append(float(els[7]))
        mag2.append(float(els[8]))
        rad.append(float(els[9]))
    return mag1, mag2, rad

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



def get_Ebv(result):
    ebvMapNorth = EBV.EBVmap()
    datadir = os.environ.get("SIMS_DUSTMAPS_DIR")
    ebvMapNorth.readMapFits(os.path.join(datadir, "DustMaps/SFD_dust_4096_ngp.fits"))
    ebvMapSouth = EBV.EBVmap()
    ebvMapSouth.readMapFits(os.path.join(datadir, "DustMaps/SFD_dust_4096_sgp.fits"))
    # Calculate galactic coordinates:
    gLon = []
    gLat = []
    for i in range(len(result['raJ2000'])):
        gcoord = afwCoord.IcrsCoord(afwGeom.Point2D(result[i]['raJ2000'], result[i]['decJ2000']), afwGeom.radians).toGalactic()
        gLon.append(gcoord.getL().asRadians())
        gLat.append(gcoord.getB().asRadians())
    return EBV.calculateEbv(gLon, gLat, ebvMapNorth, ebvMapSouth, interp = True)

def get_TotalSDSSMags(result, bandpasses=('u','g','r','i','z')):
    datadir = os.environ.get("SIMS_SED_LIBRARY_DIR")
    tpath = os.getenv('SDSS_THROUGHPUTS')
    bands = {"u":None, "g":None, "r":None, "i":None, "z":None}
    for k in bands:
        bands[k] = Bandpass()
        bands[k].readThroughput(os.path.join(tpath, "sdss_%s.dat"%k))
    # Set up phi, the wavelength-normalized system response for each filter,
    # for each bandpass for manyMagCalc method.
    bplist = []
    for f in ['u','g','r','i','z']:
        bands[f].sbTophi()
        bplist.append(bands[f])
    ids = result['galid']
    diskfile = result['sedFilenameDisk']
    bulgefile = result['sedFilenameBulge']
    agnfile = result['sedFilenameAgn']

    diskmn = result['magNormDisk']
    bulgemn = result['magNormBulge']
    agnmn = result['magNormAgn']

    bulgeAv = result['internalAvBulge']
    diskAv = result['internalAvDisk']

    redshift = result['redshift']

    imsimband = Bandpass()
    imsimband.imsimBandpass()
    sedDict = {}
    retMags = dict([(k, []) for k in bands])
    a_int = None
    b_int = None
    tmpwavelen = None
    for id, df, dm, dav, bf, bm, bav, af, am, z in zip(ids, diskfile, diskmn, diskAv, 
            bulgefile, bulgemn, bulgeAv, agnfile, agnmn, redshift):
        tmpflux = None
        for comp in ((df, dm, dav, 'galaxySED', False), (bf, bm, bav, 'galaxySED', False), (af, am, None, 'agnSED', True)):
        #Zero out the AGN contribution
        #for comp in ((df, dm, dav, 'galaxySED', False), (bf, bm, bav, 'galaxySED', False), (af, 99.99, None, 'agnSED', True)):
            if not comp[0] == u'None':
                if sedDict.has_key(comp[0]):
                    sed = copy.deepcopy(sedDict[comp[0]])
                else:
                    sed = Sed()
                    print os.path.join(datadir,comp[3],comp[0])
                    sed.readSED_flambda(os.path.join(datadir,comp[3],comp[0]))
		    if comp[4]:
		        sed.resampleSED(wavelen_match=tmpwavelen)
                    sedDict[comp[0]] = sed
                if a_int is None:
                    phiarray, dlambda = sed.setupPhiArray(bplist)
                    a_int, b_int = sed.setupCCMab()
		    #Careful, this assumes that a disk or bulge sed is read
		    #before any agn sed
		    tmpwavelen = sed.wavelen
                fNorm = sed.calcFluxNorm(comp[1], imsimband)
                sed.multiplyFluxNorm(fNorm)
                #I guess this assumes rv=3.1??
                if comp[2]:
                    sed.addCCMDust(a_int, b_int, A_v=comp[2])
		wavelenArr=sed.wavelen
		if tmpflux is None:
		    tmpflux = sed.flambda
		else:
	            tmpflux += sed.flambda
	newgal = Sed(wavelen=wavelenArr, flambda=tmpflux)
        #a_mw, b_mw = sed.setupCCMab()
        #sed.addCCMDust(a_mw, b_mw, A_v=mwav)
        newgal.redshiftSED(z, dimming=True)
	newgal.resampleSED(wavelen_match=bplist[0].wavelen)
	newgal.flambdaTofnu()
        mags = newgal.manyMagCalc(phiarray, dlambda)
        for i,k in enumerate(['u','g','r','i','z']):
            retMags[k].append(mags[i])
    return retMags

def get_TotalMags(result, bandpasses=('u','g','r','i','z','y')):
    datadir = os.environ.get("SIMS_SED_LIBRARY_DIR")
    tpath = os.getenv('LSST_THROUGHPUTS_DEFAULT')
    bands = {"u":None, "g":None, "r":None, "i":None, "z":None, "y":None}
    for k in bands:
        bands[k] = Bandpass()
        bands[k].readThroughput(os.path.join(tpath, "total_%s.dat"%k))
    # Set up phi, the wavelength-normalized system response for each filter,
    # for each bandpass for manyMagCalc method.
    bplist = []
    for f in ['u','g','r','i','z','y']:
        bands[f].sbTophi()
        bplist.append(bands[f])
    ids = result['galid']
    diskfile = result['sedFilenameDisk']
    bulgefile = result['sedFilenameBulge']
    agnfile = result['sedFilenameAgn']

    diskmn = result['magNormDisk']
    bulgemn = result['magNormBulge']
    agnmn = result['magNormAgn']

    bulgeAv = result['internalAvBulge']
    diskAv = result['internalAvDisk']

    redshift = result['redshift']

    imsimband = Bandpass()
    imsimband.imsimBandpass()
    sedDict = {}
    retMags = dict([(k, []) for k in bands])
    a_int = None
    b_int = None
    tmpwavelen = None
    for id, df, dm, dav, bf, bm, bav, af, am, z in zip(ids, diskfile, diskmn, diskAv, 
            bulgefile, bulgemn, bulgeAv, agnfile, agnmn, redshift):
        tmpflux = None
        for comp in ((df, dm, dav, 'galaxySED', False), (bf, bm, bav, 'galaxySED', False), (af, am, None, 'agnSED', True)):
        #Zero out the AGN contribution
        #for comp in ((df, dm, dav, 'galaxySED', False), (bf, bm, bav, 'galaxySED', False), (af, 99.99, None, 'agnSED', True)):
            if not comp[0] == u'None':
                if sedDict.has_key(comp[0]):
                    sed = copy.deepcopy(sedDict[comp[0]])
                else:
                    sed = Sed()
                    print os.path.join(datadir,comp[3],comp[0])
                    sed.readSED_flambda(os.path.join(datadir,comp[3],comp[0]))
		    if comp[4]:
		        sed.resampleSED(wavelen_match=tmpwavelen)
                    sedDict[comp[0]] = sed
                if a_int is None:
                    phiarray, dlambda = sed.setupPhiArray(bplist)
                    a_int, b_int = sed.setupCCMab()
		    #Careful, this assumes that a disk or bulge sed is read
		    #before any agn sed
		    tmpwavelen = sed.wavelen
                fNorm = sed.calcFluxNorm(comp[1], imsimband)
                sed.multiplyFluxNorm(fNorm)
                #I guess this assumes rv=3.1??
                if comp[2]:
                    sed.addCCMDust(a_int, b_int, A_v=comp[2])
		wavelenArr=sed.wavelen
		if tmpflux is None:
		    tmpflux = sed.flambda
		else:
	            tmpflux += sed.flambda
	newgal = Sed(wavelen=wavelenArr, flambda=tmpflux)
        #a_mw, b_mw = sed.setupCCMab()
        #sed.addCCMDust(a_mw, b_mw, A_v=mwav)
        newgal.redshiftSED(z, dimming=True)
	newgal.resampleSED(wavelen_match=bplist[0].wavelen)
	newgal.flambdaTofnu()
        mags = newgal.manyMagCalc(phiarray, dlambda)
        for i,k in enumerate(['u','g','r','i','z','y']):
            retMags[k].append(mags[i])
    return retMags

def plotHist(result, dfac, bfac, xlab, minimum=0, maximum=4.1):
    idxs = numpy.where(result['majorAxisBulge'] > 0)
    nonzerobulge = numpy.degrees(numpy.sqrt(result[idxs[0]]['majorAxisBulge']*result[idxs[0]]['minorAxisBulge']))
    idxs = numpy.where(result['majorAxisDisk'] > 0)
    nonzerodisk = numpy.degrees(numpy.sqrt(result[idxs[0]]['majorAxisDisk']*result[idxs[0]]['minorAxisDisk']))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(nonzerodisk*dfac, bins=numpy.arange(minimum, maximum, 0.1), alpha=0.5, normed=1, label='Disk')
    ax.hist(nonzerobulge*bfac, bins=numpy.arange(minimum, maximum, 0.1), alpha=0.5, normed=1, label='Bulge')
    plt.xlabel("$%s$ (arcsec)"%(xlab))
    plt.ylabel("Normalized counts")
    plt.legend()
    plt.show()

def plotEllipHist(result):
    idxs = numpy.where(result['majorAxisBulge'] > 0)
    bulge_pa = result[idxs[0]]['positionAngleBulge']
    bulge_a = numpy.degrees(result[idxs[0]]['majorAxisBulge'])
    bulge_b = numpy.degrees(result[idxs[0]]['minorAxisBulge'])
    idxs = numpy.where(result['majorAxisDisk'] > 0)
    disk_a = numpy.degrees(result[idxs[0]]['majorAxisDisk'])
    disk_b = numpy.degrees(result[idxs[0]]['minorAxisDisk'])
    disk_pa = result[idxs[0]]['positionAngleDisk']
    E = (bulge_a**2-bulge_b**2)/(bulge_a**2+bulge_b**2)
    #E = (bulge_a-bulge_b)/(bulge_a+bulge_b)
    bulge_e1 = E*numpy.cos(2*bulge_pa)
    bulge_e2 = E*numpy.sin(2*bulge_pa)
    #E = (disk_a**2-disk_b**2)/(disk_a**2+disk_b**2)
    E = (disk_a-disk_b)/(disk_a+disk_b)
    disk_e1 = E*numpy.cos(2*disk_pa)
    disk_e2 = E*numpy.sin(2*disk_pa)
    e1 = numpy.append(bulge_e1, disk_e1)
    e2 = numpy.append(bulge_e2, disk_e2)
    fig = plt.figure(figsize=(6,9))
    ax = fig.add_subplot(311)
    ax.hist(disk_e1, bins=numpy.arange(-1, 1, 0.04), alpha=0.5, normed=1, label='e1 -- Disk')
    ax.hist(disk_e2, bins=numpy.arange(-1, 1, 0.04), alpha=0.5, normed=1, label='e2 -- Disk')
    plt.ylabel("Normalized counts")
    plt.legend()
    ax = fig.add_subplot(312)
    ax.hist(bulge_e1, bins=numpy.arange(-1, 1, 0.04), alpha=0.5, normed=1, label='e1 -- Bulge')
    ax.hist(bulge_e2, bins=numpy.arange(-1, 1, 0.04), alpha=0.5, normed=1, label='e2 -- Bulge')
    plt.ylabel("Normalized counts")
    plt.legend()
    ax = fig.add_subplot(313)
    ax.hist(e1, bins=numpy.arange(-1, 1, 0.04), alpha=0.5, normed=1, label='e1 -- Combined')
    ax.hist(e2, bins=numpy.arange(-1, 1, 0.04), alpha=0.5, normed=1, label='e2 -- Combined')
    plt.ylabel("Normalized counts")
    plt.legend()
    print "done plotting... Time to show"
    plt.show()


def intNofz(z_o, z):
    return -z_o*numpy.exp(-z/z_o)*(2.*z_o**2 + 2.*z_o*z + z**2)

def coilNofz(z, z_o, norm=1., zlimits=(0.,30.)):
    area = intNofz(z_o, zlimits[1]) - intNofz(z_o, zlimits[0])
    return norm*z**2*numpy.exp(-z/z_o)/area

def get_z_o(imag):
    return 0.05*imag - 0.84

def plotNofz(rarr, iarr, redshifts, area, compfunc=coilNofz, 
        magranges_I=[(18,19), (19,20), (20,21), (21,22), (22,23), (23,24)]):
    #from coil (2004)
    z_o_corr = [0.091, 0.136, 0.197, 0.239, 0.288, 0.324]
    z_o_err = [0.01, 0.007, 0.005, 0.005, 0.006, 0.01]
    #from Lupton: http://www.sdss.org/dr5/algorithms/sdssUBVRITransform.html#Lupton2005
    derived_I = (rarr - 1.2444*(rarr-iarr) - 0.3820) #pm 0.0078
    fig = plt.figure(figsize=(15,15))
    i = 1
    magranges_I = numpy.asarray(magranges_I)
    midranges = numpy.asarray([(el[0]+el[1])/2. for el in magranges_I])
    for r,z_o,z_o_err in zip(magranges_I, z_o_corr, z_o_err):
        print z_o, midranges[i-1]
        idxs = numpy.where((derived_I > r[0])&(derived_I < r[1]))[0]
	ax=fig.add_subplot(3,2,i)
        hist = plt.hist(redshifts[idxs], bins=numpy.arange(0,3.,0.1), weights=[1./area for el in redshifts[idxs]], label="Base Catalog")
	norm = numpy.sum(hist[0])*(hist[1][1] - hist[1][0])
        cpts = numpy.arange(0, 3., 0.05)	
	plt.plot(cpts, compfunc(cpts, z_o+3*z_o_err, norm=norm), label="z$_o$+=3$\sigma$", ls="--", c='r', lw=2)
	plt.plot(cpts, compfunc(cpts, z_o, norm=norm), label="z$_o$ = z$_o$", c='k', lw=2)
	plt.plot(cpts, compfunc(cpts, z_o-3*z_o_err, norm=norm), label="z$_o$ -= 3$\sigma$", ls="-.", c='r', lw=2)
	plt.ylabel("N/deg/$3$x$10^4$s$^{-1}$km")
	plt.xlabel("redshift")
	plt.title("%i < i < %i"%(r[0], r[1]))
        if i == 1:
            plt.legend(loc="upper right")
	i += 1
    plt.tight_layout()
    plt.show()

def plotNofzResid(rarr, iarr, redshifts, compfunc=coilNofz):
    #from coil (2004)
    magranges_I = [(18,19), (19,20), (20,21), (21,22), (22,23), (23,24)]
    z_o_corr = [0.091, 0.136, 0.197, 0.239, 0.288, 0.324]
    z_o_err = [0.01, 0.007, 0.005, 0.005, 0.006, 0.01]
    #from Lupton: http://www.sdss.org/dr5/algorithms/sdssUBVRITransform.html#Lupton2005
    derived_I = (rarr - 1.2444*(rarr-iarr) - 0.3820) #pm 0.0078
    fig = plt.figure(figsize=(14,14))
    i = 1
    for r,z_o,z_o_err in zip(magranges_I, z_o_corr, z_o_err):
        idxs = numpy.where((derived_I > r[0])&(derived_I < r[1]))[0]
	ax=fig.add_subplot(3,2,i)
	hist = numpy.histogram(redshifts[idxs], bins=numpy.arange(0,3.,0.1))
	norm = numpy.sum(hist[0])*(hist[1][1] - hist[1][0])
        cpts = numpy.arange(0.05, 3., 0.05)	
	pvals = compfunc(numpy.arange(0.05, 2.9, 0.1), z_o, norm=norm)
	plt.errorbar(numpy.arange(0.05, 2.9, 0.1), hist[0]/pvals, yerr=numpy.sqrt(hist[0])/pvals)
	plt.plot(cpts, compfunc(cpts, z_o, norm=norm)/compfunc(cpts, z_o+3*z_o_err, norm=norm))
	plt.plot(cpts, compfunc(cpts, z_o, norm=norm)/compfunc(cpts, z_o, norm=norm))
	plt.plot(cpts, compfunc(cpts, z_o, norm=norm)/compfunc(cpts, z_o-3*z_o_err, norm=norm))
	plt.ylabel("Fraction")
	plt.xlabel("redshift")
	plt.title("%i < I < %i"%(r[0], r[1]))
	i += 1
    plt.show()

def plotNofzResidCum(rarr, iarr, redshifts, compfunc=coilNofz):
    #from coil (2004)
    magranges_I = [(18,19), (19,20), (20,21), (21,22), (22,23), (23,24)]
    z_o_corr = [0.091, 0.136, 0.197, 0.239, 0.288, 0.324]
    z_o_err = [0.01, 0.007, 0.005, 0.005, 0.006, 0.01]
    #from Lupton: http://www.sdss.org/dr5/algorithms/sdssUBVRITransform.html#Lupton2005
    derived_I = (rarr - 1.2444*(rarr-iarr) - 0.3820) #pm 0.0078
    fig = plt.figure(figsize=(14,14))
    i = 1
    magranges_I = numpy.asarray(magranges_I)
    midranges = numpy.asarray([(el[0]+el[1])/2. for el in magranges_I])
    for r,z_o,z_o_err in zip(magranges_I,z_o_corr, z_o_err):
        idxs = numpy.where((derived_I > r[0])&(derived_I < r[1]))[0]
	ax=fig.add_subplot(3,2,i)
	hist = numpy.histogram(redshifts[idxs], bins=numpy.arange(0,3.,0.1))
	histcum = numpy.cumsum(hist[0])
	norm = numpy.sum(hist[0])*(hist[1][1] - hist[1][0])
        cpts = numpy.arange(0.05, 3., 0.05)	
	pvals = numpy.cumsum(compfunc(numpy.arange(0.05, 2.9, 0.1), z_o, norm=norm))
	plt.errorbar(numpy.arange(0.05, 2.9, 0.1), histcum/pvals, yerr=numpy.sqrt(histcum)/pvals)
	tvals = numpy.cumsum(compfunc(cpts, z_o, norm=norm))
	plt.plot(cpts, numpy.cumsum(compfunc(cpts, z_o+z_o_err, norm=norm))/tvals)
	plt.plot(cpts, tvals/tvals)
	plt.plot(cpts, numpy.cumsum(compfunc(cpts, z_o-z_o_err, norm=norm))/tvals)
	plt.ylabel("Fraction")
	plt.xlabel("redshift")
	plt.ylim((0.8, 1.2))
	plt.title("%i < I < %i"%(r[0], r[1]))
	i += 1
    plt.show()

def loadFile(file, colTypeMap, baseCatalog=None):
    fh = open(file)
    hdr = fh.readline()
    colnames = hdr.rstrip().split(",")
    if baseCatalog:
        for name in colnames:
            if name in baseCatalog:
                continue
            else:
                raise ValueError("input base catalog does not have key %s"%(name))
    else:
        baseCatalog = {}
        for name in colnames:
            baseCatalog[name] = []
    for l in fh:
        for name, typeMap, el in zip(colnames, colTypeMap, l.rstrip().split(",")):
            baseCatalog[name].append(typeMap(el))
    fh.close()
    return baseCatalog
