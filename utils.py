from __future__ import with_statement
from collections import OrderedDict
import warnings
import numpy
from StringIO import StringIO
from sqlalchemy import types as satypes, Column, Table, Index
import string, random
#from http://stackoverflow.com/questions/2257441/python-random-string-generation-with-upper-case-letters-and-digits
def id_generator(size=8, chars=string.ascii_lowercase):
    return ''.join(random.choice(chars) for x in range(size))

def guessDtype(dataPath, numGuess, delimiter, **kwargs):
    cnt = 0
    teststr = ''
    with open(dataPath) as fh:
        while cnt < numGuess:
            teststr += fh.readline()
            cnt += 1
    dataArr = numpy.genfromtxt(StringIO(teststr), dtype=None, names=True, delimiter=delimiter, **kwargs)
    return dataArr.dtype

def buildTypeMap():
    npTypeMap = {}
    for name in dir(numpy):
        obj = getattr(numpy, name)
        if hasattr(obj, 'dtype'):
            try:
                npn = obj(0)
                nat = numpy.asscalar(npn)
                npTypeMap[npn.dtype.char] = type(nat)
            except:
                pass
    return npTypeMap

def createSQLTable(dtype, tableid, idCol, metadata, customTypeMap={}):
    npTypeMap = buildTypeMap()
    nativeTypes = OrderedDict()
    for name in dtype.names:
        try:
            nativeTypes[name] = npTypeMap[dtype[name].char]
        except KeyError:
            warings.warn("No mapping available for column %s (%s).  It will not be included in the autoloaded database.\n"\
                         "You may be able to fix this by passing a custom dtype to the class constructor."%(name, dtype[name]))
            continue
    sqlColumns = []
    for name in nativeTypes:
        try:
            sqlType = satypes._type_map[nativeTypes[name]]
        except KeyError:
            warnings.warn("Could not find a mapping to a SQL type for native type: %s"%(nativeTypes[name]))
        #Look for string like types
        if name in customTypeMap:
            sqlColumns.append(Column(name, customTypeMap[name], primary_key=(idCol==name)))
        elif 'format' in dir(nativeTypes[name]):
            sqlColumns.append(Column(name, type(sqlType)(dtype[name].itemsize), primary_key=(idCol==name)))
        else:
            sqlColumns.append(Column(name, sqlType, primary_key=(idCol==name)))
    if tableid is None:
        tableid = id_generator()
    datatable = Table(tableid, metadata, *sqlColumns)
    metadata.create_all()
    return datatable

def loadTable(datapath, datatable, delimiter, dtype, engine, indexCols=[], skipLines=1, chunkSize=100000, **kwargs):
    cnt = 0
    with open(datapath) as fh:
        while cnt < skipLines:
            fh.readline()
            cnt += 1
        cnt = 0
        tmpstr = ''
        for l in fh:
            tmpstr += l
            cnt += 1
            if cnt%chunkSize == 0:
                print "Loading chunk #%i"%(int(cnt/chunkSize))
                dataArr = numpy.genfromtxt(StringIO(tmpstr), dtype=dtype, delimiter=delimiter, **kwargs)
                engine.execute(datatable.insert(), [dict((name, numpy.asscalar(l[name])) for name in l.dtype.names) for l in dataArr])
                tmpstr = ''
        #Clean up the last chunk
        if len(tmpstr) > 0:
            dataArr = numpy.genfromtxt(StringIO(tmpstr), dtype=dtype, delimiter=delimiter, **kwargs)
            engine.execute(datatable.insert(), [dict((name, numpy.asscalar(l[name])) for name in l.dtype.names) for l in dataArr])

    for col in indexCols:
    	if hasattr(col, "__iter__"):
        	print "Creating index on %s"%(",".join(col))
        	colArr = (datatable.c[c] for c in col)
        	i = Index('%sidx'%''.join(col), *colArr)
        else:
        	print "Creating index on %s"%(col)
        	i = Index('%sidx'%col, datatable.c[col])

        i.create(engine)

def loadData(dataPath, dtype, delimiter, tableId, idCol, engine, metaData, numGuess, **kwargs):
	if dtype is None:
		dtype = guessDtype(dataPath, numGuess, delimiter)
	dataTable = createSQLTable(dtype, tableId, idCol, metaData)
	loadTable(dataPath, dataTable, delimiter, dtype, engine, **kwargs)
	return dataTable.name
