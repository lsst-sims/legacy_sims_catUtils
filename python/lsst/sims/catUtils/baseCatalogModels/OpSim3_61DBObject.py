from lsst.sims.catalogs.db import DBObject
from lsst.sims.utils import ObservationMetaData
from .BaseCatalogModels import BaseCatalogConfig
from lsst.utils import getPackageDir
from collections import OrderedDict
import numpy
import math
import os

from lsst.sims.utils import CircleBounds

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
    raColName = 'fieldradeg'
    decColName = 'fielddecdeg'

    # columnNames maps column names as stored in the OpSim3_61 database
    # to column names as they occur in the OpSim v3.3.5 summary table
    columnNames = [('obsHistID','obshistid'),
                   ('SIM_SEED','expdate'),
                   ('pointingRA','fieldra'),
                   ('pointingDec','fielddec'),
                   ('moonRA','moonra'),
                   ('moonDec','moondec'),
                   ('rotSkyPos','rotskypos'),
                   ('rotTelPos','rottelpos'),
                   ('rawSeeing','rawseeing'),
                   ('sunAlt','sunalt'),
                   ('moonAlt','moonalt'),
                   ('dist2Moon','dist2moon'),
                   ('moonPhase','moonphase'),
                   ('expMJD','expmjd'),
                   ('visitExpTime','exptime'),
                   ('airmass', 'airmass'),
                   ('filter', 'filter')]

    #columnTypes is a dict indicating the data types
    #of the columns referred to in columnNames.
    #columns not mentioned are floats
    columnTypes ={'SIM_SEED':int,
                  'filter':(str,1),
                  'obsHistID':numpy.int64}

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

    def getObservationMetaData(self, searchCenter, searchRadius, constraint=None,
                               fovRadius=1.75, makeCircBounds=True, makeBoxBounds=False):
        """
        This method will perform a query on the OpSim3_61 for pointings within a given region
        on the sky.  It will then return a list of ObservationMetaData corresponding to the
        OpSim pointings returned by that query.

        param [in] searchCenter is a tuple (RA, Dec) indicating the center of the region
        of the sky wherein to search for OpSim pointings (in degrees).

        param [in] searchRadius is the radius of the region of the sky wherein to search
        for OpSim pointings (in degrees).

        param [in] constraint an optional constraint on another column.  Should be a string
        of the form 'airmass=1.1'

        param [in] fovRadius is the radius (in degrees) of the resulting ObservationMetaDatas'
        fields of view.  If you have specified a box bounds this will be half the length of each
        side of the box.

        param [in] makeCircBounds a boolean telling this method to construct an ObservationMetaData
        for a circular region of the sky centered on the specified OpSim pointing

        param [in] make BoxBounds a boolean telling this method to construct an
        ObservationMetaData for a square region of the sky centered on the specified
        OpSim pointing
        """

        spatialBound = CircleBounds(numpy.radians(searchCenter[0]), numpy.radians(searchCenter[1]),
                                    numpy.radians(searchRadius))

        # get the results of the query as a numpy recarray
        result = self.executeConstrainedQuery(spatialBound, constraint=constraint)

        obs_list = []
        ra_list = numpy.degrees(result['pointingRA'])
        dec_list = numpy.degrees(result['pointingDec'])
        rotSkyPos_list = numpy.degrees(result['rotSkyPos'])
        for pointing, ra, dec, rotSkyPos in zip(result, ra_list, dec_list, rotSkyPos_list):

            mjd = pointing['expMJD']

            #because makeCircBounds defaults to True, check whether the user is
            #requesting boxBounds before deciding to instantiate
            #circBounds

            boundType = None
            boundLength = None
            if makeBoxBounds:
                boundType = 'box'
                boundLength = numpy.array([fovRadius/math.cos(numpy.radians(dec)), fovRadius])
            elif makeCircBounds:
                boundType = 'circle'
                boundLength = fovRadius
            else:
                raise ValueErr("Need either makeBoxBounds or makeCircBounds")

            obs = ObservationMetaData(pointingRA=ra, pointingDec=dec,
                                      rotSkyPos=rotSkyPos, mjd=mjd,
                                      bandpassName=result['filter'][0],
                                      boundType=boundType, boundLength=boundLength)

            raw_opsim_dict = dict([(k, pointing[k]) for k in result.dtype.names])
            obs.OpsimMetaData = raw_opsim_dict
            obs_list.append(obs)

        return obs_list
