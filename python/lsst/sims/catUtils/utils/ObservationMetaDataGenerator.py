import os
import eups
from lsst.sims.catalogs.generation.db import DBObject

class ObservationMetaDataGenerator(object):

    def __init__(self, address=None):
        if address is None:
            dbPath = os.path.join(eups.productDir('sims_data'),'OpSimData')
            self.address=os.path.join('sqlite://',dbPath,'opsimblitz1_1133_sqlite.db')
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
                             ('airmass','airmass',float,'airmass')]

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
                               visitExpTime=None, airmass=None)

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

        """

        query = self.baseQuery+ ' FROM SUMMARY WHERE'

        for element in self.columnMapping:
            value = eval(element[3])
            if value is not None:
                if isinstance(value,tuple):
                    if len(value)>2:
                        raise RuntimeError('Cannot pass a tuple longer than 2 elements to getObservationMetaData: %s is len %d' \
                                           % (element[3], len(value)))

                    query += ' %s between %s and %s' % (element[0], value[0], value[1])
                else:
                    query += ' %s = %s' % (element[0], value)

        self.opsimdb.execute_arbitrary(query, dtype=self.dtype)
