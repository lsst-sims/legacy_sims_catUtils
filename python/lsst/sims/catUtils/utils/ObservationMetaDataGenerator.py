import os
import eups
import numpy
from collections import OrderedDict
from lsst.sims.catalogs.generation.db import DBObject, ObservationMetaData

__all__ = ["ObservationMetaDataGenerator"]

class ObservationMetaDataGenerator(object):
    
    def _put_quotations(self, val):
        if val[0]=='\'':
            return val
        else:
            return "'%s'" % val

    def __init__(self, address=None):
        if address is None:
            dbPath = os.path.join(eups.productDir('sims_data'),'OpSimData/')
            self.address='sqlite:///' + dbPath + 'opsimblitz1_1133_sqlite.db'
        else:
            self.address=address

        self.opsimdb = DBObject(address=self.address)

        #27 January 2015
        #self.columnMapping is an OrderedDict of tuples.  The keys of the dict are how users will
        #refer to the OpSim summary table columns (ie. how they are called in getObservationMetaDAta).
        #The 0th element of the tuple is how the column is named in the OpSim db.
        #The 1st element of the tuple is how PhoSim refers to the quantity (as of OpSim3_61DBObject.py)
        #The 2nd element of the tuple is the datatype of the column
        #The 3d element of the tuple is any coordinate transformation required between user interface
        #and the OpSim database (i.e. OpSim stores all angles in radians; we would like users to be
        #able to specify angles in degrees)
        #
        #Note that this conforms to an older
        #PhoSim API.  At some time in the future, both this and OpSim3_61DBObject.py
        #(and possibly phoSimCatalogExamples.py) will need to be updated to
        #reflect what PhoSim actually expects now.
        self.columnMapping = OrderedDict([('obsHistID',('obsHistID','Opsim_obshistid',numpy.int64,None)),
                                         ('expDate',('expDate','SIM_SEED',int, None)),
                                         ('fieldRA',('fieldRA','Unrefracted_RA',float,numpy.radians)),
                                         ('fieldDec',('fieldDec','Unrefracted_Dec',float,numpy.radians)),
                                         ('moonRA',('moonRA','Opsim_moonra',float, numpy.radians)),
                                         ('moonDec',('moonDec','Opsim_moondec',float,numpy.radians)),
                                         ('rotSkyPos',('rotSkyPos','Opsim_rotskypos',float,numpy.radians)),
                                         ('telescopeFilter',('filter','Opsim_filter',(str,1),self._put_quotations)),
                                         ('rawSeeing',('rawSeeing','Opsim_rawseeing',float,None)),
                                         ('sunAlt',('sunAlt','Opsim_sunalt',float, numpy.radians)),
                                         ('moonAlt',('moonAlt','Opsim_moonalt',float,numpy.radians)),
                                         ('dist2Moon',('dist2Moon','Opsim_dist2moon',float, numpy.radians)),
                                         ('moonPhase',('moonPhase','Opsim_moonphase',float,None)),
                                         ('expMJD',('expMJD','Opsim_expmjd',float,None)),
                                         ('altitude',('altitude','Opsim_altitude',float,numpy.radians)),
                                         ('azimuth',('azimuth','Opsim_azimuth',float,numpy.radians)),
                                         ('visitExpTime',('visitExpTime','exptime',float,None)),
                                         ('airmass',('airmass','airmass',float,None)),
                                         ('m5',('fiveSigmaDepth',None,float,None)),
                                         ('skyBrightness',('filtSkyBrightness',None,float,None))])

        dtypeList = []
        self.baseQuery = 'SELECT'
        for column in self.columnMapping:
            dtypeList.append((self.columnMapping[column][0],self.columnMapping[column][2]))
            if self.baseQuery != 'SELECT':
                self.baseQuery += ','
            self.baseQuery += ' ' + self.columnMapping[column][0]

        self.dtype = numpy.dtype(dtypeList)


    def getObservationMetaData(self, obsHistID=None, expDate=None, fieldRA=None, fieldDec=None,
                               moonRA=None, moonDec=None, rotSkyPos=None, telescopeFilter=None,
                               rawSeeing=None, sunAlt=None, moonAlt=None, dist2Moon=None,
                               moonPhase=None, expMJD=None, altitude=None, azimuth=None,
                               visitExpTime=None, airmass=None, skyBrightness=None,
                               m5=None, boundType='circle', boundLength=0.1, limit=None):

        """
        This method will query the OpSim database according to user-specified constraints
        and return a list of of ObservationMetaData instantiations consistent with those
        constraints.

        The parameters that can be passed in are all either tuples of the form (min,max) or
        straight values corresponding to the quantities in the OpSim database that can be
        contrained.  If tuples are passed, the query will ask for pointings such that
        min < x < max for the quantity being contrained.  If a number is passed, then
        the query will ask that x == value.

        The quantities that can be constrained are:

        fieldRA in degrees
        fieldDec in degrees
        altitude in degrees
        azimuth in degrees

        moonRA in degrees
        moonDec in degrees
        moonAlt in degrees
        moonPhase (a value from 1 to 100 indicating how much of the moon is illuminated)
        dist2Moon the distance between the telescope pointing and the moon in degrees

        sunAlt in degrees

        rotSkyPos (the angle of the sky with respect to the camera coordinate system) in degrees
        telescopeFilter (u,g,r,i,z,y)

        airmass
        rawSeeing

        visitExpTime the exposure time in seconds
        obsHistID the integer used by OpSim to label pointings
        expDate is the date of the exposure (units????)
        expMJD is the MJD of the exposure
        m5 is the five sigma depth of the observation
        skyBrightness

        You must also specify the size of the fields of view returned in the
        ObservationMetaData instantiations by specifying boundType (default 'circle')
        and boundLength (in degrees; default a radius of 0.1).   See documentation in
        sims_catalogs_generation/../db/spatialBounds.py for more details.

        limit limts the number of rows returned

        """

        query = self.baseQuery+ ' FROM SUMMARY'
        
        nConstraints = 0 #the number of constraints in this query

        for column in self.columnMapping:
            value = eval(column)
            if value is not None:
                if nConstraints > 0:
                    query += ' AND'
                else:
                    query += ' WHERE '
                    
                if isinstance(value,tuple):
                    if len(value)>2:
                        raise RuntimeError('Cannot pass a tuple longer than 2 elements '+
                                           'to getObservationMetaData: %s is len %d'
                                           % (column, len(value)))

                    if self.columnMapping[column][3] is not None:
                        vmin = self.columnMapping[column][3](value[0])
                        vmax = self.columnMapping[column][3](value[1])
                    else:
                        vmin = value[0]
                        vmax = value[1]

                    query += ' %s > %s AND %s < %s' % \
                             (self.columnMapping[column][0], vmin, self.columnMapping[column][0], vmax)
                else:
                    if self.columnMapping[column][3] is not None:
                        vv = self.columnMapping[column][3](value)
                    else:
                        vv = value
                    query += ' %s == %s' % (self.columnMapping[column][0], vv)
                
                nConstraints += 1
        
        if limit is not None:
            query += ' LIMIT %d' % limit

        if nConstraints==0 and limit is None:
            raise RuntimeError('You did not specify any contraints on your query;' +
                               ' you will just return ObservationMetaData for all poitnings')

        results = self.opsimdb.execute_arbitrary(query, dtype=self.dtype)
        
        obs_output = []
        
        for pointing in results:
            phoSimMetadata=OrderedDict()
            for column in self.columnMapping:
                if self.columnMapping[column][1] is not None:
                    phoSimMetadata[self.columnMapping[column][1]] = (pointing[self.columnMapping[column][0]],
                                                                     pointing[self.columnMapping[column][0]].dtype)
            
            obs_metadata = ObservationMetaData(m5=pointing['fiveSigmaDepth'],
                                               boundType=boundType, boundLength=boundLength,
                                               phoSimMetadata=phoSimMetadata)
            
            obs_output.append(obs_metadata)
        
        return obs_output
        
