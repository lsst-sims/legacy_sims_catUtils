"""
This script will illustrate how to write your own InstanceCatalog daughter
class.  Recall that the InstanceCatalog class is defined in

sims_catalogs_measures/python/lsst/sims/catalogs/measures/instance/InstanceCatalog.py

After reading this script, you should look at the classes defined in
the scripts in

sims_catUtils/python/lsst/sims/catUtils/exampleCatalogDefinitions/

These files contain most of the physically interesint InstanceCatalog
daughter classes defined thus far.  This will show you how the functionality
described below is implemented in actual use cases.

Because InstanceCatalog has its own __init__() method that needs to be called
to set up member variables on which the catalog relies (and because we do not
want to be wed to python's Super() API), InstanceCatalog daughter classes should
not define their own __init__() methods.  Any member variables unique to the
daughter class should be set in the class definition outside of any method.

The most important member variable of any InstanceCatalog daughter class is
the list column_outputs.  column_outputs is a list of the names of the columns
that are to be output to the catalog.

InstanceCatalog and its daughter classes find the data to fill in these columns in
three ways.  In this order, they:

1) assess whether or not they have a method to calculate the column value

2) assess whether or not they can directly query the column from their
associated CatalogDBObject

3) assess whether they know any default values for those columns.

(1) and (3) must be defined in the InstanceCatalog daughter class.
(2) is provided by the CatalogDBObject that has been passed to __init__() when
the InstanceCatalog daughter class is instantiated.

(1) relies on getters.  Literally, the InstanceCatalog daughter class inspects
its methods and looks for a method named get_columnName (whatever columnName
actually is).  These methods can perform whatever calculations are necessary
to calculate the values of their returned columns.  They must return numpy
arrays of the calculated values.

Columns calculated by getters can rely on the values of other columns.  These
are accessed using the method self.column_by_name('nameOfColumn'), which returns
a numpy array of column values.  i.e.

def get_myNewColumn(self):
    old = self.column_by_name('myOldColumn')
    ....some math that depends on the values of myOldColumn...
    return numpy.array([valuesOfNewColumn])

It is possible to calculate multiple columns in a single getter using the
@compound decorator defined in

sims_catalogs_measures/python/lsst/sims/catalogs/measures/instance/decorators.py

We we will see an example of this decorator below.  One calls it with the names
of the columns to be calculated, i.e.

@compound('col1', 'col2')
def get_manyColumns(self):
    ...some math...
    return numpy.array([col1Values, col2Values])

In this case, the getter must return an 2-D numpy array in which each row is a different
set of column values and the column values occur in the same order that they occur
in @compound()

(3) -- setting column values from a default -- relies on the member variable
default_columns.  default_columns is a list of tuples.  Each tuple consists of
three elements: the name of the column being defaulted; the default value of
that column; the datatype of that column.

The script below defines the InstanceCatalog daughter class
TutorialCatalog.  This catalog contains getters for the sum, difference
and quotient of raJ2000 and decJ2000 (which we assume are supplied by the
CatalogDBObject).  It sets default values for the dummy columns df1 and df2.

At the end, this script writes an instantiation of TutorialCatalog to
tutorial_catalog.txt
"""

import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
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
obs_metadata = ObservationMetaData(unrefractedRA=220.0, unrefractedDec=19.0,
                                   boundType='circle', boundLength=0.1)

cat = TutorialCatalog(myDB, obs_metadata=obs_metadata)
cat.write_catalog('tutorial_catalog.txt')
