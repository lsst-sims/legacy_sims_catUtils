"""
This is a version of tutorial04_CatalogDBObject.py without the running commentary
"""

import os
import numpy
import sqlite3

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
myDB = TutorialDB(driver='sqlite', database='tutorialDB.db')



print 'First show all of the columns in the raw database'
myDB.show_db_columns()
print '\n'




print 'Then show all of the columns in the CatalogDBObject'
myDB.show_mapped_columns()
print '\n'

print 'now do a rough, by-hand query of the columns (this returns all of the rows)'

#specify the names of the columns to be queried
colNames = ['rowNumber', 'TwiceColumn1', 'col1', 'ThreeTimesColumn2', 'col2']

result = myDB.query_columns(colnames=colNames, chunk_size=5)

for chunk in result:
    for row in chunk:
        print row
print '\n'


print 'now apply a constraint to the query'
result = myDB.query_columns(colnames=colNames, constraint='id<6', chunk_size=5)
for chunk in result:
    for row in chunk:
        print row
print '\n'

print 'now we will write a cartoon catalog class to write out this cartoon database'

from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached

class TestCatalog(InstanceCatalog):
    column_outputs = ['col1', 'col2', 'col3', 'col4', 'threeMinusFour']
    catalog_type = 'test_catalog'

    @cached
    def get_threeMinusFour(self):
        three = self.column_by_name('col3')
        four = self.column_by_name('col4')
        return three-four

myCat = TestCatalog(myDB)
myCat.write_catalog('cartoon_catalog.txt')
