import numpy
#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catalogs.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils import ValidationUtils, SersicUtils
import argparse

def galaxyCountValidation(dbobj, doSaveFigs, fileBase):
    binsize = 0.5
    bins = numpy.arange(14.1, 28.6, binsize)
    compx, compy, comperr = ValidationUtils.read_Durham_counts("baseValidation/test_data/idata.txt")
    compx = numpy.array(compx)
    compy = numpy.array(compy)
    comperr = numpy.array(comperr)
    compx += 0.6
    comp_means = numpy.histogram(compx, bins=bins, weights=compy)[0]/numpy.histogram(compx, bins=bins)[0]
    norm = (2.*2.25)*(2.*2.25)
    ii = 0
    basemean = []
    binx = []
    for i in bins[:-1]:
        qstr = "select count(*) from %s where i_ab_orig between %f and %f"%(dbobj.tableid, i, i+binsize)
        cur = dbobj.engine.execute(qstr)
        result = cur.fetchone()
        basemean.append(result[0]/norm)
        binx.append((i+i+binsize)/2.)
        print i, result[0], result[0]/(result[0]/norm/comp_means[ii]) - result[0], result[0]/norm, comp_means[ii], result[0]/norm/comp_means[ii]
        ii += 1
    if doSaveFigs:
        ValidationUtils.plotN(compx, compy, comperr, basemean, binx, filename=fileBase+"_galaxyCounts.png")
        ValidationUtils.plotNResidCume(compx, compy, comperr, basemean, binx, filename=fileBase+"_galaxyCountsResidCume.png")
    else:
        ValidationUtils.plotN(compx, compy, comperr, basemean, binx)
        ValidationUtils.plotNResidCume(compx, compy, comperr, basemean, binx)
    
def redshiftDistributionValidation(dbobj, obs_metadata, doSaveFigs, fileBase):
    dmag = 1.
    magbins = numpy.arange(18., 25., dmag)
    dz = 0.1
    zbins = numpy.arange(0., 3.+dz, dz)
    zbinx = [(el1+el2)/2. for el1,el2 in zip(zbins[:-1], zbins[1:])]
    #These are the values for z_o from Coil (2004) for the mag ranges above
    z_o_corr = [0.091, 0.136, 0.197, 0.239, 0.288, 0.324]
    z_o_err = [0.01, 0.007, 0.005, 0.005, 0.006, 0.01]

    #z_o_corr = 0.0417*numpy.array([(el1+el2)/2. for el1,el2 in zip(magbins[:-1], magbins[1:])]) - 0.744
    histData = []
    area = (2.*2.25)*(2.*2.25)
    for i in magbins:
        #from Lupton: http://www.sdss.org/dr5/algorithms/sdssUBVRITransform.html#Lupton2005
        t = dbobj.getCatalog("ref_catalog_galaxy", 
                             obs_metadata=obs_metadata, 
                             constraint="i_ab-0.1 between %f and %f"%(i, i+dmag))
                             #constraint="r_ab - 1.24444*(r_ab-i_ab) - 0.3820 between %f and %f"%(i,i+dmag))
        hist = None
        ii = 0
        for chunk, colMap in t.iter_catalog_chunks(chunk_size=100000):
            redshift = chunk[colMap['redshift']]
            thist = numpy.histogram(redshift, bins=zbins, normed=False)
            print "Done %i -- %f to %f"%(ii, i, i+dmag )
            if hist is None:
                hist = thist[0]
            else:
                hist += thist[0]
            ii += 1
        histData.append(hist)
    if doSaveFigs:
        ValidationUtils.plotNofz(histData, zbinx, magbins, z_o_corr, z_o_err, area, filename=fileBase+"_galaxyNofz.png")
    else:
        ValidationUtils.plotNofz(histData, zbinx, magbins, z_o_corr, z_o_err, area)

def ellipticityDistributionValidation(dbobj, obs_metadata, doSaveFigs, fileBase):
    t = dbobj.getCatalog("ref_catalog_galaxy",
                         obs_metadata=obs_metadata,
                         constraint="i_ab < 25.5 and e1 is NULL")
    hist = None
    diskSersic = SersicUtils(1.)
    bulgeSersic = SersicUtils(4.)
    ntot = 10000
    e1s = []
    e2s = []
    r_hls = []
    fh = open("galaxy_ellip.dat", "w")
    fh.write("id, e1, e2, Ixx, Iyy, Ixy, r1, r_hl\n")
    for chunk, colMap in t.iter_catalog_chunks(chunk_size=10):
        diskmag = chunk[colMap['DiskLSSTi']]
        diskmajor = chunk[colMap['majorAxisDisk']]
        diskminor = chunk[colMap['minorAxisDisk']]
        diskpa = chunk[colMap['positionAngleDisk']]

        bulgemag = chunk[colMap['BulgeLSSTi']]
        bulgemajor = chunk[colMap['majorAxisBulge']]
        bulgeminor = chunk[colMap['minorAxisBulge']]
        bulgepa = chunk[colMap['positionAngleBulge']]
        galid = chunk[colMap['galid']]
        data = numpy.array(zip(galid, diskmag, numpy.radians(diskmajor), numpy.radians(diskminor), numpy.radians(diskpa), bulgemag, numpy.radians(bulgemajor), numpy.radians(bulgeminor), numpy.radians(bulgepa)), 
                           dtype=[('galid', (str, 50)), ('diskmag', float), ('diskmajor', float), ('diskminor', float), ('diskpa', float), 
                                  ('bulgemag', float), ('bulgemajor', float), ('bulgeminor', float), ('bulgepa', float)])
        for i in xrange(len(data)):
            nbulge, ndisk = SersicUtils.getNphots(ntot, data['bulgemag'][i], data['diskmag'][i])
            dx, dy = diskSersic.calcDistribution(data['diskmajor'][i], data['diskminor'][i], data['diskpa'][i], ndisk)
            bx, by = bulgeSersic.calcDistribution(data['bulgemajor'][i], data['bulgeminor'][i], data['bulgepa'][i], nbulge)
            e1, e2, Ixx, Iyy, Ixy, r1, r_hl = SersicUtils.getEllip(numpy.concatenate((dx, bx)), numpy.concatenate((dy, by)))
            fh.write(",".join([str(el) for el in (data['galid'][i], e1, e2, numpy.degrees(Ixx), numpy.degrees(Iyy), 
                                                  numpy.degrees(Ixy), numpy.degrees(r1), numpy.degrees(r_hl))])+"\n")
            e1s.append(e1)
            e2s.append(e2)
            r_hls.append(r_hls)
    fh.close()
    import pdb;pdb.set_trace()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--skipCounts", help="Skip running the galaxy count validation", 
                        action="store_true")
    parser.add_argument("--skipNofz", help="Skip running the redshift distribution validation", 
                        action="store_true")
    parser.add_argument("--skipEllip", help="Skip running the ellipticity distribution validation", 
                        action="store_true")
    parser.add_argument("--savePlots", help="Save plots as PNG files instead of displaying them?",
                        action="store_true")
    parser.add_argument("--baseOutputName", help="Base filename for output plots",
                        default='validation')
    parser.add_argument("--boxsize", help="size of the side of the box to get in degrees",
                        default=8., type=float)
    args = parser.parse_args()
    obs_metadata = ObservationMetaData(circ_bounds=dict(ra=0., dec=0., radius=args.boxsize/2.))
    dbobj = CatalogDBObject.from_objid('galaxyBase')
    filename = None
    if not args.skipCounts:
        galaxyCountValidation(dbobj, args.savePlots, args.baseOutputName)
    if not args.skipNofz:
        redshiftDistributionValidation(dbobj, obs_metadata, args.savePlots, args.baseOutputName)
    if not args.skipEllip:
        ellipticityDistributionValidation(dbobj, obs_metadata, args.savePlots, args.baseOutputName)
