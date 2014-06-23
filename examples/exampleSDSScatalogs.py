from lsst.sims.photUtils.photometry import PhotometryStars PhotometryGalaxies
from lsst.sims.catalogs.measures.instance import compound

class sdssStars(PhotometryStars):
    
    @compound('sdss_uRecalc', 'sdss_gRecalc', 'sdss_rRecalc', 
              'sdss_iRecalc', 'sdss_zRecalc',
              'sdss_uBulge', 'sdss_gBulge', 'sdss_rBulge', 'sdss_iBulge', 'sdss_zBulge',
              'sdss_uDisk', 'sdss_gDisk', 'sdss_rDisk', 'sdss_iDisk', 'sdss_zDisk',
              'sdss_uAgn', 'sdss_gAgn', 'sdss_rAgn', 'sdss_iAgn', 'sdss_zAgn')
    def get_all_sdss_mags(self):
        """
        example getter for sdss galaxy magnitudes
        
        bandPassRoot is the root of the names of the files in which
        the bandpasses are stored
        
        """
        idNames = self.column_by_name('galid')
        bandPassList = ['u','g','r','i','z','y']
        return self.meta_magnitudes_getter(idNames, bandPassList, bandPassRoot = 'sdss_')



class sdssGalaxies(PhotometryGalaxies):
  
    @compound('sdss_u','sdss_g','sdss_r','sdss_i','sdss_z')
    def get_sdss_magnitudes(self):
        """
        example getter for sdss stellar magnitudes
        
        bandPassRoot is the root of the names of the files in which
        the bandpasses are stored
        """
        idNames = self.column_by_name('id')
        bandPassList = ['u','g','r','i','z']
        bandPassRoot = 'sdss_'
        return self.meta_magnitudes_getter(idNames, bandPassList, bandPassRoot = bandPassRoot)
