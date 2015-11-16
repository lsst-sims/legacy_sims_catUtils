"""
This is a version of tutorial01.ipynb without the
running commentary
"""

import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached, compound

class TutorialCatalog(InstanceCatalog):

    #this the member variable that defines what columns should be output
    #when we call write_catalog
    column_outputs = ['raJ2000', 'decJ2000', 'sum', 'difference', 'quotient',
                      'df1', 'df2']

    #this is the member variable that contains default values for columns
    #that are neither defined in the database nor set by a getter
    default_columns = [('df1', 3.5, float), ('df2', 'default', (str,7))]

    @cached
    def get_quotient(self):
        """
        A getter for the single column 'quotient'

        The @cached decorator means that, whenever this getter is called,
        its values are stored to be accessed any time self.column_by_name('quotient')
        is called.  Without @cached, self.column_by_name('quotient') would have
        to call get_quotient()
        """
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')

        #because self.column_by_name() returns a numpy array, we can use
        #numpy array's matrix formalism to calculate the result
        return ra/dec

    @compound('sum', 'difference')
    def get_sumAndDifference(self):
        """
        A compound getter for the columns 'sum' and 'difference'

        @compound automatically applies @cached
        """
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')

        #note that the columns must be returned in the order in which
        #they are declared to @compound()
        return numpy.array([ra+dec, ra-dec])

myDB = CatalogDBObject.from_objid('allstars')
obs_metadata = ObservationMetaData(pointingRA=220.0, pointingDec=19.0,
                                   boundType='circle', boundLength=0.1,
                                   mjd=57388.0)

cat = TutorialCatalog(myDB, obs_metadata=obs_metadata)
cat.write_catalog('tutorial_catalog.txt')
