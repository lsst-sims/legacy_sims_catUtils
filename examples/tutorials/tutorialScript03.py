"""
This tutorial will start to examine what goes into a CatalogDBObject
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
    tableid = 'example'
    idColKey = 'id'
    objid = 'tutorial_DBobject'
    
    columns = [('2xc1','2.0*col1'),
              ('3xc2','3.0*col2'),
              ('rowNumber','col5',str,10)]

makeTestDB('tutorialDB.db', size=10)
myDB = TutorialDB(address='sqlite:///tutorialDB.db')

print 'First show all of the columns in the database'
myDB.show_db_columns()
print '\n'

print 'Then show all of the columns in the CatalogDBObject'
myDB.show_mapped_columns()
print '\n'

print 'now do a rough, by-hand query of the columns (this returns all of the columns)'
colNames = ['rowNumber', '2xc1', 'col1']
result = myDB.query_columns(colnames=colNames, chunk_size=5)

for chunk in result:
    for row in chunk:
        print row
