from lsst.sims.catalogs.generation.db import ChunkIterator, CatalogDBObject, ObservationMetaData
from lsst.sims.catalogs.generation.db.spatialBounds import SpatialBounds
from collections import OrderedDict
import numpy
import math

__all__ = ["OpSim3_61DBObject"]

class OpSim3_61DBObject(CatalogDBObject):
    """Meta Data Database Object base class

    """
    objid = 'opsim3_61'

    #: This is the default address.  Simply change this in the class definition for other
    #: endpoints.
    dbAddress = "mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST"

    tableid = 'output_opsim3_61'
    objectTypeId = -1
    generateDefaultColumnMap = False
    #: Note that identical observations may have more than one unique
    #: obshistid, so this is the id, but not for unique visits.
    #: To do that, group by expdate.
    idColKey = 'Opsim_obshistid'
    bandColKey = 'Opsim_filter'
    raColKey = 'Unrefracted_RA'
    decColKey = 'Unrefracted_Dec'
    #: These are interpreted as SQL strings.
    raColName = 'fieldra*180./PI()'
    decColName = 'fielddec*180./PI()'
    mjdColName = 'expmjd'
    columns = [('SIM_SEED', 'expdate', int),
               ('Unrefracted_RA', 'fieldra'),
               ('Unrefracted_Dec', 'fielddec'),
               ('Opsim_moonra', 'moonra'),
               ('Opsim_moondec', 'moondec'),
               ('Opsim_rotskypos', 'rotskypos'),
               ('Opsim_rottelpos', 'rottelpos'),
               ('Opsim_filter', 'filter', str, 1),
               ('Opsim_rawseeing', 'rawseeing'),
               ('Opsim_sunalt', 'sunalt'),
               ('Opsim_moonalt', 'moonalt'),
               ('Opsim_dist2moon', 'dist2moon'),
               ('Opsim_moonphase', 'moonphase'),
               ('Opsim_obshistid', 'obshistid', numpy.int64),
               ('Opsim_expmjd', 'expmjd'),
               ('Opsim_altitude', 'altitude'),
               ('Opsim_azimuth', 'azimuth'),
               ('exptime', 'exptime'),
               ('airmass', 'airmass')]
               
    def __init__(self, address=None):
        self.defaultColumns = [col[0] for col in self.columns]
        super(OpSim3_61DBObject, self).__init__(address=address)

    def getObjectTypeId(self):
        raise NotImplementedError("Metadata has no object type")

    def getSpatialModel(self):
        raise NotImplementedError("Metadata has no spatial model")

    def getObservationMetaData(self, obshistid, radiusDeg, makeCircBounds=True, makeBoxBounds=False, colnames=None):
        if colnames is None:
            colnames = self.defaultColumns
        chunks = self.query_columns(colnames=colnames, constraint="obshistid=%i"%obshistid)
        #The query will only return one row (hopefully)
        result = chunks.next()
        ra = result[self.raColKey][0]
        dec = result[self.decColKey][0]
        
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
                           unrefractedRA=math.degrees(ra), unrefractedDec=math.degrees(dec),
                           boundLength=boundLength,
                           phoSimMetadata=OrderedDict([(k, (result[k][0], result[k][0].dtype)) for k in result.dtype.names]))

    def query_columns(self, colnames=None, chunk_size=None,
                      boundType=None, boundLength=None,
                      unrefractedRA=None, unrefractedDec=None,
                      mjd_bounds=None, constraint=None):
        """Execute a query

        **Parameters**

            * colnames : list or None
              a list of valid column names, corresponding to entries in the
              `columns` class attribute.  If not specified, all columns are
              queried.
            * chunk_size : int (optional)
              if specified, then return an iterator object to query the database,
              each time returning the next `chunk_size` elements.  If not
              specified, all matching results will be returned.
            * circ_bounds : dict (optional)
              a dictionary with the keys 'ra', 'dec', and 'radius' measured, in
              degrees
            * box_bounds : dict (optional)
              a dictionary with the keys 'ra_min', 'ra_max', 'dec_min', 'dec_max',
              measured in degrees
            * mjd_bounds : dict (optional)
              a dictionary with the keys 'mjd_min' and 'mjd_max' used to bound the 
              query in time.
            * constraint : str (optional)
              if constraint exists it will be used verbatim as a filter on the query 

        **Returns**

            * result : structured array or iterator
              If chunk_size is not specified, then result is a structured array whose
              columns names are specified in the 'columns' class variable. If chunk_size 
              is specified, then result is an iterator over structured arrays of the 
              given size.

        """
        query = self._get_column_query(colnames)

        bounds = None

        if boundType is not None:
            if unrefractedRA is None or unrefractedDec is None or boundLength is None:
                raise RuntimeError("in Opsim3_61DBObject query_columns, bound improperly specified")
            bounds = SpatialBounds.getSpatialBounds(boundType,unrefractedRA,unrefractedDec,boundLength)

        query = self.filter(query, bounds)
        if mjd_bounds is not None:
            query = query.filter("%s between %f and %f"%(self.mjdColName,
                                  mjd_bounds['mjd_min'], mjd_bounds['mjd_max']))
        if constraint is not None:
            query = query.filter(constraint)
        return ChunkIterator(self, query, chunk_size)
