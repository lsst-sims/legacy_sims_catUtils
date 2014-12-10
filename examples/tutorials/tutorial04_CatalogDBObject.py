"""
This tutorial will start to examine what goes into a CatalogDBObject

After reading this script, you should probably look at the classes defined
by the scripts in

sims_catUtils/python/lsst/sims/catUtils/baseCatalogModels/

These files will contain most of the CatalogDBObject daughter classes that
have so far been defined.  Here you can see how the functionality described
below is implemented for actual, physically interesting cases.
"""

import os
import numpy
import sqlite3

"""
First, we require some code that will make a database.  The  code below
will make a database containing a table named 'example' that contains
nonsense columns col1-col5"
"""

def makeTestDB(filename='tutorialDatabase.db', size=1000, seedVal=None, **kwargs):
    """
    make a test database for tutorial purposes
    @param size: Number of rows in the database
    @param seedVal: Random seed to use
    """

    if os.path.exists(filename):
        os.unlink(filename)

    conn = sqlite3.connect(filename)
    c = conn.cursor()
    try:
        c.execute('''CREATE TABLE example
                     (id int, col1 int, col2 int, col3 real, col4 real, col5 text)''')
        conn.commit()
    except:
        raise RuntimeError("Error creating database.")
    if seedVal:
        seed(seedVal)

    numpy.random.seed(seedVal)
    c1 = numpy.random.randint(0,100,size)
    c2 = numpy.random.randint(0,100,size)
    c3 = numpy.random.random_sample(size)
    c4 = numpy.random.random_sample(size)

    for i in xrange(size):
        c5 = '%sth_row' % i
        qstr = '''INSERT INTO example VALUES (%i, %i, %i, %f, %f,
                     '%s')'''%\
                   (i,c1[i],c2[i],c3[i],c4[i],c5)
        c.execute(qstr)

    conn.commit()
    conn.close()


"""
Now we will discuss how to define a daughter class of CatalogDBobject.
Recall that CatalogDBobject is defined in

sims_catalogs_generation/python/lsst/sims/catalogs/generation/db/dbConnnection.py

As with InstanceCatalog, daughter classes of CatalogDBobject should not contain
an __init__().  Instead, any member variables that are unique to the daughter
class should be defined in the class definition.

The most important member variables in a CatalogDBObject are

dbAddress -- this is the connection string to the database.  This can be passed
             in through __init__()

tableid -- this is a string indicating with what table in the database this
           CatalogDBObject is associated

idColKey -- this is a string indicating the name of the column that serves as
            a unique identifier for objects in the database

objid -- like catalog_type in InstanceCatalog() daughter classes, this is the
         handle that is passed to CatalogDBObject.from_objid() when you want
         to automatically instantiate this CatalogDBObject daughter class.

The purpose of a CatalogDBObject is to provide InstanceCatalog with a connection
to the database.  Every InstanceCatalog instantiation has a member variable
self.db_obj.  This is a CatalogDBObject.  If you examine the method
write_catalog in InstanceCatalog, you will see that InstanceCatalog gets its
data by calling self.db_obj.query_columns().  This supplies all of the columns
not defined by getters or default values in the InstanceCatalog. query_columns()
can get two different kinds of columns.  It can get columns that simply exist
in the database by querying the names of columns as they are stored in the
database.  It can also get columns that are simple transformations of columns
stored in the database.  This process is handled using the CatalogDBObject
member variable columns.

columns is a list.  It is somewhat analogous to the default_columns member
variable of InstanceCatalog.  Each entry in the list is a tuple.  The first
element of each tuple is the name of the transformed column.  The second
element is an expression calculating that column from a column stored in the
database.  The (optional) third element is the datatype of the transformed
column.

Consider the class StarBase from

sims_catUtils/python/lsst/sims/catUtils/baseCatalogModels/StarModels.py

It contains

    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mura/(1000.*3600.))*PI()/180.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', unicode, 40)]

Ths list tells the Catalog DBObject that:

    'id' is just 'simobjid' cast as an int

    'raJ2000' is 'ra' converted from degrees to radians (the type is implicit)

    ditto for 'decJ2000', 'glon', and 'glat'

    'magNorm' is determined from 'flux_scale' using the relationship between
    magnitude and flux

    'properMotionRa' and 'properMotionDec' are calculated by converting
    'mura' and 'mudecl' from milliarcseconds/yr to radians/yr

    'galacticAv' is 'ebv' multiplied by 3.1 and cast as a float

    'radialVelocity' is 'vrad' (no change except the column's name)

    'variabilityParameters' is a string of maximum 256 characters and is
    identical to 'varParamStr'

    'sedFilename' is a 40 character unicode string and is identical to
    'sedfilename'

Thus, there are two types of 'columns' supplied by a CatalogDBObject:
the raw contents of the database, and the raw contents of the database
plus the mapped columns.

Fortunately, there are methods to tell you what columns are available from
each CatalogDBObject.

show_db_columns() will print the raw database columns to the screen.

show_mapped_columns() will print both the raw database columns and the mapped
columns to the screen.

We demonstrate those methods below.  After that, we show an example call to
query_columns.  There will be more documentation where that actually occurs.
"""

from lsst.sims.catalogs.generation.db import CatalogDBObject

class TutorialDB(CatalogDBObject):
    """
    Here we define a tutorial CatalogDBObject
    """

    #the database table that this class queries
    tableid = 'example'

    #the column that uniquely identifies each object in the catalog
    idColKey = 'id'

    #this is the handle that can be passed to CatalogDBObject.from_objid()
    #as in tutorial00
    objid = 'tutorial_DBobject'

    #here we convert the raw database columns into mapped columns
    columns = [('TwiceColumn1','2.0*col1'),
              ('ThreeTimesColumn2','3.0*col2'),
              ('rowNumber','col5',str,10)]




makeTestDB('tutorialDB.db', size=10)
myDB = TutorialDB(address='sqlite:///tutorialDB.db')



print 'First show all of the columns in the raw database'
myDB.show_db_columns()
print '\n'




print 'Then show all of the columns in the CatalogDBObject'
myDB.show_mapped_columns()
print '\n'


"""
Now we will use the query_columns method to query the database.
Ordinarily, the user should not have to do this.  InstanceCatalog does it
automatically.  This is just for illustrative purposes.

query_columns returns a ChunkIterator as definee in

sims_catalogs_generation/python/lsst/sims/catalogs/generation/db/dbConnection.py

which returns chunks of columns at a time (the number of rows in each chunk is
specified by the chunk_size kwarg passed to query_columns)
"""


print 'now do a rough, by-hand query of the columns (this returns all of the rows)'

#specify the names of the columns to be queried
colNames = ['rowNumber', 'TwiceColumn1', 'col1', 'ThreeTimesColumn2', 'col2']

result = myDB.query_columns(colnames=colNames, chunk_size=5)

for chunk in result:
    for row in chunk:
        print row
print '\n'


"""
One also has the option of passing a constraint to query_columns to limit the
rows that are returned by the query.  In many cases, this constraint is supplied
by the ObservationMetaData (specifically the ObservationMetaData.bounds member variable)
"""

print 'now apply a constraint to the query'
result = myDB.query_columns(colnames=colNames, constraint='id<6', chunk_size=5)
for chunk in result:
    for row in chunk:
        print row
print '\n'
