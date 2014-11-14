"""
Here we demonstrate using actual mixins that have already been defined
to calculate physically real values.

Many mixins exist through out the stack.  The two biggest collections of them
are in

sims_coordUtils/python/lsst/sims/coordUtils/Astrometry.py  -- astrometry-related code
sims_photUtils/python/lsst/sims/photUtils/Photometry.py -- photometry-related code
sims_photUtils/python/lsst/sims/photUtils/EBV.py -- dust-related code
sims_photUtils/python/lsst/sims/photUtils/Variability.py -- variability-related code

Below, we demonstrate an InstanceCatalog class that uses the photometry
and astrometry mixins to calculate physically meaningful column values
"""

import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData, haversine
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached
from lsst.sims.coordUtils import AstrometryStars
from lsst.sims.photUtils import PhotometryStars

class TutorialCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):

    #lsst_u is defined in PhotometryStars in Photometry.py; it is the magnitude in the LSST u band
    #ra/decObserved are defined in AstrometryStars i Astrometry.py; they are the RA and Dec corrected
    #               for precession, nutation, aberration, and refraction
    #shift is a column defined by a getter below; it is the difference between (raJ2000, decJ2000) and
    #             (raObserved, decObserved)
    column_outputs = ['raJ2000', 'decJ2000', 'lsst_u', 'raObserved', 'decObserved',
                      'shift']

    #transform all of the angles into degrees
    transformations = {'raJ2000':numpy.degrees, 'decJ2000':numpy.degrees,
                      'raObserved':numpy.degrees, 'decObserved':numpy.degrees,
                       'shift':numpy.degrees}

    #a handle to be passed to CatalogDBObject.getCatalog() (see tutorial02)
    catalog_type = 'tutorial_catalog'

    @cached
    def get_shift(self):
        """
        A getter for the angular distance between the unrefracted raJ2000, decJ2000
        and the corrected raObserved, decObserved

        Note that because all angles are handled inside of the stack as radians,
        the returned angular distance will also be in radians
        """
        r0 = self.column_by_name('raJ2000')
        d0 = self.column_by_name('decJ2000')
        r1 = self.column_by_name('raObserved')
        d1 = self.column_by_name('decObserved')
        
        return haversine(r0, d0, r1, d1)




#write the catalog directly
myDB = CatalogDBObject.from_objid('allstars')
obs_metadata = ObservationMetaData(unrefractedRA=220.0, unrefractedDec=19.0,
                                   boundType='circle', boundLength=0.1, 
                                   mjd=52000.0)

cat = TutorialCatalog(myDB, obs_metadata=obs_metadata)
cat.write_catalog('tutorial_astrometry_photometry.txt')



#write the catalog using CatalogDBObject.getCatalog()
obs_metadata = ObservationMetaData(unrefractedRA=120.0, unrefractedDec=-5.0,
                                   boundType='circle', boundLength=0.1, 
                                   mjd=52000.0)

cat = myDB.getCatalog('tutorial_catalog', obs_metadata=obs_metadata)
cat.write_catalog('tutorial_astrometry_photometry_get.txt')
