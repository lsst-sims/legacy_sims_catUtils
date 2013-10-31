from dbConnection import ChunkIterator, DBObject
from collections import OrderedDict
import numpy
import math

class ObservationMetaData(object):
    """Observation Metadata
    
    This class contains any metadata for a query which is associated with
    a particular telescope pointing, including bounds in RA and DEC, and
    the time of the observation.

    **Parameters**

        * circ_bounds : dict (optional)
          a dictionary with the keys 'ra', 'dec', and 'radius' measured, in
          degrees
        * box_bounds : dict (optional)
          a dictionary with the keys 'ra_min', 'ra_max', 'dec_min', 'dec_max',
          measured in degrees
        * mjd : float (optional)
          The MJD of the observation..
        * metadata : dict (optional)
          a dictionary containing arbitrary metadata

    **Examples**::

        >>> data = ObservationMetaData(dict(('ra_min', 0), ('ra_max', 10), ('dec_min', 10), ('dec_max', 20)))

    """
            
    def __init__(self, circ_bounds=None, box_bounds=None, mjd=None, metadata=None):
        if circ_bounds is not None and box_bounds is not None:
            raise ValueError("Passing both circ_bounds and box_bounds")
        self.circ_bounds = circ_bounds
        self.box_bounds = box_bounds
        self.mjd = mjd
        self.metadata = metadata

class MetaDataDBObject(DBObject):
    """Meta Data Database Object base class

    """
    objid = 'opsim3_61'
    tableid = 'output_opsim3_61'
    objectTypeId = -1
    spatialModel = 'none'
    #: Note that identical observations may have more than one unique
    #: obshistid, so this is the id, but not for unique visits.
    #: To do that, group by expdate.
    idColKey = 'Opsim_obshistid'
    bandColKey = 'Opsim_filter'
    raColKey = 'Unrefracted_RA'
    decColKey = 'Unrefracted_Dec'
    #: These are interpreted as SQL strings.
    raColName = 'fieldra'
    decColName = 'fielddec'
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
               ('Opsim_azimuth', 'azimuth')]
    def __init__(self, address=None):
        self.defaultColumns = [col[0] for col in self.columns]
        super(MetaDataDBObject, self).__init__(address=address)

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
        if makeCircBounds:
            circ_bounds = dict(ra=math.degrees(ra), 
                                    dec=math.degrees(dec), 
                                    radius=radiusDeg)
            box_bounds = None
        elif makeBoxBounds:
            box_bounds = dict(ramin=math.degrees(ra)-radiusDeg/math.cos(dec), 
                              ramax=math.degrees(ra)+radiusDeg/math.cos(dec),
                              decmin=math.degrees(dec)-radiusDeg,
                              decmax=math.degrees(dec)+radiusDeg)
            circ_bounds = None
        else:
            raise ValueErr("Need either circ_bounds or box_bounds")
        return ObservationMetaData(circ_bounds=circ_bounds, box_bounds=box_bounds, 
                           metadata=OrderedDict([(k, (result[k][0], result[k][0].dtype)) for k in result.dtype.names]))

    def query_columns(self, colnames=None, chunk_size=None,
                      circ_bounds=None, box_bounds=None, mjd_bounds=None,
                      constraint=None):
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

        query = self.filter(query, circ_bounds=circ_bounds, 
                    box_bounds=box_bounds)
        if mjd_bounds is not None:
            query = query.filter("%s between %f and %f)"%(self.mjdColName,
                                  mjd_bounds['mjd_min'], mjd_bounds['mjd_max']))
        if constraint is not None:
            query = query.filter(constraint)
        return ChunkIterator(self, query, chunk_size)
