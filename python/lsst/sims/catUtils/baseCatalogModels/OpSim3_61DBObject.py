from lsst.sims.catalogs.db import DBObject
from lsst.sims.utils import ObservationMetaData
from .BaseCatalogModels import BaseCatalogConfig
from lsst.utils import getPackageDir
from collections import OrderedDict
import numpy
import math
import os

__all__ = ["OpSim3_61DBObject"]


class OpSim3_61DBObject(DBObject):
    """
    This is a class which allows you to query the databse of OpSim runs (v 3_61)
    """
    config = BaseCatalogConfig()
    config.load(os.path.join(getPackageDir("sims_catUtils"), "config", "db.py"))
    host = config.host
    port = config.port
    database = config.database
    driver = config.driver

    #This is the table that this DBObject is querying
    tableid = 'output_opsim3_61'

    #: These are interpreted as SQL strings.
    #They are passed to SpatialBounds objects
    #and indicate the RA and Dec columns in degrees
    raColName = 'fieldra*180./PI()'
    decColName = 'fielddec*180./PI()'

    #columnNames maps column names as stored in the database
    #to column names as are required by PhoSim.
    #For each tuple, the first string is the name
    #required by PhoSim.  The second string is
    #the name stored in the database
    columnNames = [('Opsim_obshistid','obshistid'),
               ('SIM_SEED','expdate'),
               ('pointingRA','fieldra'),
               ('pointingDec','fielddec'),
               ('Opsim_moonra','moonra'),
               ('Opsim_moondec','moondec'),
               ('Opsim_rotskypos','rotskypos'),
               ('Opsim_rottelpos','rottelpos'),
               ('Opsim_filter','filter'),
               ('Opsim_rawseeing','rawseeing'),
               ('Opsim_sunalt','sunalt'),
               ('Opsim_moonalt','moonalt'),
               ('Opsim_dist2moon','dist2moon'),
               ('Opsim_moonphase','moonphase'),
               ('Opsim_expmjd','expmjd'),
               ('Opsim_altitude','altitude'),
               ('Opsim_azimuth','azimuth'),
               ('exptime','exptime'),
               ('airmass','airmass')]

    #columnTypes is a dict indicating the data types
    #of the columns referred to in columnNames.
    #columns not mentioned are floats
    columnTypes ={'SIM_SEED':int,
                  'Opsim_filter':(str,1),
                  'Opsim_obshistid':numpy.int64}

    def __init__(self, driver=None, host=None, port=None, database=None):
        super(OpSim3_61DBObject, self).__init__(driver=driver, host=host, port=port, database=database)

    def getBasicQuery(self):
        """
        This class creates a string that is an SQL query of all the
        columns needed to create an ObservationMetaData suitable for
        generating a PhoSim input catalog.  No constraints
        are placed on the query (it will grab all entries in the database)
        """
        query = 'SELECT '
        for name in self.columnNames:
            if query != 'SELECT ':
                query += ', '
            query += name[1]
        query += ' FROM ' + self.tableid

        dtype=numpy.dtype([(name[0], self.columnTypes[name[0]]) if name[0] in self.columnTypes else (name[0], float) for name in self.columnNames])

        return query, dtype

    def executeConstrainedQuery(self, spatialBound, constraint=None):
        """
        This method will perform a query which is contrained on the sky by spatialBound.
        You can add a supplementary constraint on another column using the 'contraint' kwarg.

        param [in] spatialBound a SpatialBound object specifying the region of the sky from
        which you want your pointings drawn

        param [in] constraint an optional constraint on another column.  Should be a string
        of the form 'airmass=1.1'

        param [out] a recarray of all of the pointings matching your constraints.  This recarray
        will contain all of the columns necessary for constructing an ObservationMetaData object
        suitable for generating a PhoSim input catalog.
        """
        query, dtype = self.getBasicQuery()
        query += ' WHERE '+spatialBound.to_SQL(self.raColName, self.decColName)
        if constraint is not None:
            query += ' and %s' % constraint


        results = self.execute_arbitrary(query, dtype=dtype)
        return results

    def getObservationMetaData(self, obshistid, radiusDeg, makeCircBounds=True, makeBoxBounds=False):
        """
        This class creates an ObservationMetaData object suitable for generating a PhoSim input
        catalog based on an actual OpSim pointing.

        param [in] obshistid is the obshistid of the pointing to be used

        param [in] radiusDeg is the radius (in degrees) of the desired region on the sky.
        If you have specified a box bounds this will be half the length of each side of the box.

        param [in] makeCircBounds a boolean telling this method to construct an ObservationMetaData
        for a circular region of the sky centered on the specified OpSim pointing

        param [in] make BoxBounds a boolean telling this method to construct an
        ObservationMetaData for a square region of the sky centered on the specified
        OpSim pointing
        """
        query, dtype = self.getBasicQuery()
        query += ' WHERE obshistid = %i' % obshistid

        result = self.execute_arbitrary(query, dtype=dtype)

        ra = result['pointingRA'][0]
        dec = result['pointingDec'][0]

        #because makeCircBounds defaults to True, check whether the user is
        #requesting boxBounds before deciding to instantiate
        #circBounds

        boundType = None
        boundLength = None
        if makeBoxBounds:
            boundType = 'box'
            boundLength = numpy.array([radiusDeg/math.cos(dec),radiusDeg])
        elif makeCircBounds:
            boundType = 'circle'
            boundLength = radiusDeg
        else:
            raise ValueErr("Need either makeBoxBounds or makeCircBounds")

        return ObservationMetaData(boundType=boundType,
                                   boundLength=boundLength,
                                   phoSimMetaData=OrderedDict([(k, (result[k][0], result[k][0].dtype)) for k in result.dtype.names]))
