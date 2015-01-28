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
        #columnMapping is a list of tuples.  The first entry in each tuple is what the column
        #are called by OpSim.  The second value in each tuple is what OpSim3_61DBObject.py
        #calls the same quantity.  The third value is its datatype.  The fourth value is
        #how the column will be referred to in the args of getObservationMetaData (this is there
        #because 'filter' means something special to python, so we have to use 'telescopeFilter
        #in the args of getObservationMetaData.
        #
        #Note that this conforms to an older
        #PhoSim API.  At some time in the future, both this and OpSim3_61DBObject.py
        #(and possibly phoSimCatalogExamples.py) will need to be updated to
        #reflect what PhoSim actually expects now.
        self.columnMapping = [('obsHistID','Opsim_obshistid',numpy.int64, 'obsHistID'),
                             ('expDate','SIM_SEED',int, 'expDate'),
                             ('fieldRA','Unrefracted_RA',float, 'fieldRA'),
                             ('fieldDec','Unrefracted_Dec',float, 'fieldDec'),
                             ('moonRA','Opsim_moonra',float, 'moonRA'),
                             ('moonDec','Opsim_moondec',float, 'moonDec'),
                             ('rotSkyPos','Opsim_rotskypos',float, 'rotSkyPos'),
                             ('filter','Opsim_filter',(str,1), 'telescopeFilter'),
                             ('rawSeeing','Opsim_rawseeing',float, 'rawSeeing'),
                             ('sunAlt','Opsim_sunalt',float, 'sunAlt'),
                             ('moonAlt','Opsim_moonalt',float, 'moonAlt'),
                             ('dist2Moon','Opsim_dist2moon',float, 'dist2Moon'),
                             ('moonPhase','Opsim_moonphase',float, 'moonPhase'),
                             ('expMJD','Opsim_expmjd',float, 'expMJD'),
                             ('altitude','Opsim_altitude',float, 'altitude'),
                             ('azimuth','Opsim_azimuth',float, 'azimuth'),
                             ('visitExpTime','exptime',float, 'visitExpTime'),
                             ('airmass','airmass',float,'airmass'),
                             ('fiveSigmaDepth','m5',float,'m5'),
                             ('filtSkyBrightness','skyBrightness',float,'skyBrightness')]

        self.columnUnitTransformations = {'fieldRA':numpy.radians, 'fieldDec':numpy.radians,
                                          'moonRA':numpy.radians, 'moonDec':numpy.radians,
                                          'rotSkyPos':numpy.radians, 'sunAlt':numpy.radians,
                                          'moonAlt':numpy.radians, 'dist2Moon':numpy.radians,
                                          'altitude':numpy.radians, 'azimuth':numpy.radians,
                                          'filter':self._put_quotations}

        dtypeList = []
        self.baseQuery = 'SELECT'
        for element in self.columnMapping:
            dtypeList.append((element[0],element[2]))
            if self.baseQuery != 'SELECT':
                self.baseQuery += ','
            self.baseQuery += ' ' + element[0]

        self.dtype = numpy.dtype(dtypeList)


    def getObservationMetaData(self, obsHistID=None, expDate=None, fieldRA=None, fieldDec=None,
                               moonRA=None, moonDec=None, rotSkyPos=None, telescopeFilter=None,
                               rawSeeing=None, sunAlt=None, moonAlt=None, dist2Moon=None,
                               moonPhase=None, expMJD=None, altitude=None, azimuth=None,
                               visitExpTime=None, airmass=None, skyBrightness=None,
                               m5=None, boundType='circle', boundLength=0.1):

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
        moonPhase (a value from 0 to 1 indicating how much of the moon is illuminated)
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

        """

        query = self.baseQuery+ ' FROM SUMMARY WHERE'
        nWhereClauses = 0

        for element in self.columnMapping:
            value = eval(element[3])
            if value is not None:
                if nWhereClauses > 0:
                    query += ' and'
                if isinstance(value,tuple):
                    if len(value)>2:
                        raise RuntimeError('Cannot pass a tuple longer than 2 elements to getObservationMetaData: %s is len %d' \
                                           % (element[3], len(value)))

                    if element[0] in self.columnUnitTransformations:
                        vmin = self.columnUnitTransformations[element[0]](value[0])
                        vmax = self.columnUnitTransformations[element[0]](value[1])
                    else:
                        vmin = value[0]
                        vmax = value[1]

                    query += ' %s > %s and %s < %s' % (element[0], vmin, element[0], vmax)
                else:
                    if element[0] in self.columnUnitTransformations:
                        vv = self.columnUnitTransformations[element[0]](value)
                    else:
                        vv = value
                    query += ' %s == %s' % (element[0], vv)
                
                nWhereClauses += 1

        results = self.opsimdb.execute_arbitrary(query, dtype=self.dtype)
        
        obs_output = []
        
        for pointing in results:
            phoSimMetadata=OrderedDict()
            for column in self.columnMapping:
                phoSimMetadata[column[1]] = (pointing[column[0]], column[2])
        
            obs_metadata = ObservationMetaData(m5=pointing['fiveSigmaDepth'],
                                               boundType=boundType, boundLength=boundLength,
                                               phoSimMetadata=phoSimMetadata)
            
            obs_output.append(obs_metadata)
        
        return obs_output
        
