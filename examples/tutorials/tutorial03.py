"""
This is a version of tutorial03.ipynb without the running commentary
"""

import numpy
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils import haversine
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import cached
from lsst.sims.catUtils.mixins import AstrometryStars, PhotometryStars

class TutorialCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):

    column_outputs = ['raJ2000', 'decJ2000', 'lsst_u', 'raObserved', 'decObserved',
                      'shift']
    #to see where lsst_u comes from, see the PhotometryStars class in
    #sims_catUtils/python/lsst/sims/catUtils/mixins/PhotometryMixin.py
    #
    #to see where raObserved and decObserved come from, see the AstrometryStars class in
    #sims_catUtils/python/lsst/sims/catUtils/mixins/AstrometryMixin.py


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
obs_metadata = ObservationMetaData(pointingRA=220.0, pointingDec=19.0,
                                   boundType='circle', boundLength=0.1,
                                   mjd=52000.0)

cat = TutorialCatalog(myDB, obs_metadata=obs_metadata)
cat.write_catalog('tutorial_astrometry_photometry.txt')



#write the catalog using CatalogDBObject.getCatalog()
obs_metadata = ObservationMetaData(pointingRA=120.0, pointingDec=-5.0,
                                   boundType='circle', boundLength=0.1,
                                   mjd=52000.0)

cat = myDB.getCatalog('tutorial_catalog', obs_metadata=obs_metadata)
cat.write_catalog('tutorial_astrometry_photometry_get.txt')
