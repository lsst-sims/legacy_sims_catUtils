from lsst.sims.catalogs.generation.db import ChunkIterator, DBObject, ObservationMetaData
from lsst.sims.catalogs.generation.db.spatialBounds import SpatialBounds
from collections import OrderedDict
import numpy
import math

__all__ = ["OpSim3_61DBObject"]

class OpSim3_61DBObject(DBObject):
    """Meta Data Database Object base class

    """
    #: This is the default address.  Simply change this in the class definition for other
    #: endpoints.
    dbAddress = "mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST"

    tableid = 'output_opsim3_61'
    
    #: These are interpreted as SQL strings.
    raColName = 'fieldra*180./PI()'
    decColName = 'fielddec*180./PI()'

    columnNames = [('Opsim_obshistid','obshistid'),
               ('SIM_SEED','expdate'),
               ('Unrefracted_RA','fieldra'),
               ('Unrefracted_Dec','fielddec'),
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
    
    columnTypes ={'SIM_SEED':int,
                  'Opsim_filter':(str,1),
                  'Opsim_obshistid':numpy.int64}
               
    def __init__(self, address=None):
        super(OpSim3_61DBObject, self).__init__(address=address)
    
    def getBasicQuery(self):
        query = 'SELECT '
        for name in self.columnNames:
            if query != 'SELECT ':
                query += ', '
            query += name[1]
        query += ' FROM ' + self.tableid
        
        dtype=numpy.dtype([(name[0], self.columnTypes[name[0]]) if name[0] in self.columnTypes else (name[0], float) for name in self.columnNames])

        return query, dtype
   
    def executeConstrainedQuery(self, spatialBound, constraint=None):
        query, dtype = self.getBasicQuery()
        query += ' WHERE '+spatialBound.to_SQL(self.raColName, self.decColName)
        if constraint is not None:
            query += ' and %s' % constraint
       
        
        results = self.execute_arbitrary(query, dtype=dtype)
        return results
        
    def getObservationMetaData(self, obshistid, radiusDeg, makeCircBounds=True, makeBoxBounds=False, colnames=None):
        
        query, dtype = self.getBasicQuery()
        query += ' WHERE obshistid = %i' % obshistid
        
        result = self.execute_arbitrary(query, dtype=dtype)
        
        ra = result['Unrefracted_RA'][0]
        dec = result['Unrefracted_Dec'][0]
        
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
