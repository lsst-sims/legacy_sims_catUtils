import os
import numpy
import lsst.utils
from lsst.sims.catalogs.generation.db import DBObject
from lsst.sims.utils import ObservationMetaData

__all__ = ["ObservationMetaDataGenerator"]


class ObservationMetaDataGenerator(object):
    """
    A class that allows the user to generate instantiations of
    `lsst.sims.utils.ObservationMetaData` corresponding to OpSim pointings.
    The functionality includes:
    - getOpSimRecords : obtain OpSim records matching the intersection of user
        specified ranges on each column in the OpSim output database. The
        records are in the form of a `numpy.recarray`
    - ObservationMetaDataFromPointing : convert an OpSim record for a single
        OpSim Pointing to an instance of ObservationMetaData usable by catsim
        and PhoSim Instance Catalogs.
    - getObservationMetaData : Obtain a list of ObservationMetaData instances
        corresponding to OpSim pointings matching the intersection of user
        specified ranges on each column in the OpSim output database.

    The major method is ObservationMetaDataGenerator.getObservationMetaData()
    which accepts bounds on columns of the opsim summary table and returns
    a list of ObservationMetaData instantiations that fall within those
    bounds.
    """

    def __init__(self, database, driver='sqlite', host=None, port=None):
        """
        Constructor for the class

        Parameters
        ----------
        database : string, mandatory
            absolute path to the output of the OpSim database
        driver : string, optional, defaults to 'sqlite'
            driver/dialect for the SQL database
        host : hostName, optional, defaults to None,
            hostName, None is good for a local database
        port : hostName, optional, defaults to None,
            port, None is good for a local database

        Returns
        ------
        Instance of the ObserverMetaDataGenerator class

        ..notes : For testing purposes a small OpSim database is available at
        `os.path.join(getPackageDir('sims_data'), 'OpSimData/opsimblitz1_1133_sqlite.db')`
        """
        self.driver = driver
        self.host = host
        self.port = port
        self.database = database
        self.opsimdb = DBObject(driver=self.driver, database=self.database,
                                host=self.host, port=self.port)

        # 27 January 2016
        # Detect whether the OpSim db you are connecting to uses 'finSeeing'
        # as its seeing column (deprecated), or FWHMeff, which is the modern
        # standard
        summary_columns = self.opsimdb.get_column_names('Summary')

        if 'FWHMeff' in summary_columns:
            self._seeing_column = 'FWHMeff'
        else:
            self._seeing_column = 'finSeeing'

        # 27 January 2015
        # self.columnMapping is an list of tuples.  Each tuple corresponds to a column in the opsim
        # database's summary table.

        # The 0th element of the tuple is how users will
        # refer to the OpSim summary table columns (ie. how they are called in getObservationMetaDAta).
        #
        # The 1st element of the tuple is how the column is named in the OpSim db.
        #
        # The 2nd element of the tuple is how PhoSim refers to the quantity (as of OpSim3_61DBObject.py)
        # (this is None if PhoSim does not expect the quantity)
        #
        # The 3rd element of the tuple is the datatype of the column
        #
        # The 4th element of the tuple is any coordinate transformation required to go from the user interface
        # to the OpSim database (i.e. OpSim stores all angles in radians; we would like users to be
        # able to specify angles in degrees)
        #
        # The 5th column is any coordinate transformation required to go from the OpSim database
        # to the PhoSim header
        #
        # Note that this conforms to an older
        # PhoSim API.  At some time in the future, both this and OpSim3_61DBObject.py
        # (and possibly phoSimCatalogExamples.py) will need to be updated to
        # reflect what PhoSim actually expects now.
        #

        # Just use the static method to assign this attribute
        # to the mapping of names of phoSim to OpSim columns
        self.columnMapping = self.OpSimColumnMap(self._seeing_column)

        #Set up self.dtype containg the dtype of the recarray we expect back from the SQL query.
        #Also setup baseQuery which is just the SELECT clause of the SQL query
        #
        #self.active_columns will be a list containing the subset of columnMapping columns
        #that actually exist in this opsim database
        dtypeList = []
        self.baseQuery = 'SELECT'
        self.active_columns = []
        for column in self.columnMapping:
            if column[1] in summary_columns:
                self.active_columns.append(column)
                dtypeList.append((column[1],column[3]))
                if self.baseQuery != 'SELECT':
                    self.baseQuery += ','
                self.baseQuery += ' ' + column[1]

        self.dtype = numpy.dtype(dtypeList)

    @staticmethod
    def OpSimColumnMap(seeing_column):
        """
        Return a list of tuples.  Each tuple has the following elements :
        (phosimName, OpSimName, unknown, dtype, callable)

        Needs some work:
        """
        v = [('obsHistID', 'obsHistID', 'Opsim_obshistid', numpy.int64, None, None),
             ('expDate', 'expDate', 'SIM_SEED', int, None, None),
             ('fieldRA', 'fieldRA', 'Unrefracted_RA', float, numpy.radians, numpy.degrees),
             ('fieldDec', 'fieldDec', 'Unrefracted_Dec', float, numpy.radians, numpy.degrees),
             ('moonRA', 'moonRA', 'Opsim_moonra', float, numpy.radians, numpy.degrees),
             ('moonDec', 'moonDec', 'Opsim_moondec', float, numpy.radians, numpy.degrees),
             ('rotSkyPos', 'rotSkyPos', 'Opsim_rotskypos', float, numpy.radians, numpy.degrees),
             ('telescopeFilter', 'filter', 'Opsim_filter', (str, 1), lambda x: '\'{}\''.format(x), None),
             ('rawSeeing', 'rawSeeing', 'Opsim_rawseeing', float, None, None),
             ('seeing', seeing_column, None, float, None, None),
             ('sunAlt', 'sunAlt', 'Opsim_sunalt', float, numpy.radians, numpy.degrees),
             ('moonAlt', 'moonAlt', 'Opsim_moonalt', float, numpy.radians, numpy.degrees),
             ('dist2Moon', 'dist2Moon', 'Opsim_dist2moon', float, numpy.radians, numpy.degrees),
             ('moonPhase', 'moonPhase', 'Opsim_moonphase', float, None, None),
             ('expMJD', 'expMJD', 'Opsim_expmjd', float, None, None),
             ('altitude', 'altitude', 'Opsim_altitude', float, numpy.radians, numpy.degrees),
             ('azimuth', 'azimuth', 'Opsim_azimuth', float, numpy.radians, numpy.degrees),
             ('visitExpTime', 'visitExpTime', 'exptime', float, None, None),
             ('airmass', 'airmass', 'airmass', float, None, None),
             ('m5', 'fiveSigmaDepth', None, float, None, None),
             ('skyBrightness', 'filtSkyBrightness', None, float, None, None)]
        return v

    def getOpSimRecords(self, obsHistID=None, expDate=None, fieldRA=None,
                        fieldDec=None, moonRA=None, moonDec=None,
                        rotSkyPos=None, telescopeFilter=None, rawSeeing=None,
                        seeing=None, sunAlt=None, moonAlt=None, dist2Moon=None,
                        moonPhase=None, expMJD=None, altitude=None,
                        azimuth=None, visitExpTime=None, airmass=None,
                        skyBrightness=None, m5=None, boundType='circle',
                        boundLength=0.1, limit=None):
        """
        This method will query the summary table in the `self.opsimdb` database
        according to constraints specified in the input ranges and return a
        `numpy.recarray` containing the records that match those constraints. If limit
        is used, the first N records will be returned in the list.

        Parameters
        ----------
        obsHistID, expDate, fieldRA, fieldDec, moonRa, moonDec, rotSkyPos,
        telescopeFilter, rawSeeing, seeing, sunAlt, moonAlt, dist2Moon,
        moonPhase, expMJD, altitude, azimuth, visitExpTime, airmass,
        skyBrightness, m5 : tuples of length 2, optional, defaults to None
            each of these variables represent a single column (perhaps through
            an alias) in the OpSim database, and potentially in a different unit.
            if not None, the variable self.columnMapping is used to constrain
            the corresponding column in the OpSim database to the ranges specified
            in the tuples, after a unit transformation if necessary.

            The ranges must be specified in the tuple in degrees for all angles in this
            (moonRa, moonDec, rotSkyPos, sunAlt, moonAlt, dist2Moon, altitude,
            azimuth). The times in  (expMJD, are in units of MJD). visitExpTime has
            units of seconds since the start of the survey. moonPhase is a number
            from 0., to 100.
        boundType : `sims.utils.ObservationMetaData.boundType`, optional, defaults to 'circle'
            {'circle', 'box'} denoting the shape of the pointing. Further
            documentation `sims.catalogs.generation.db.spatialBounds.py``
        boundLength : float, optional, defaults to 0.1
            sets `sims.utils.ObservationMetaData.boundLenght`
        limit : integer, optional, defaults to None
            if not None, denotes max number of records returned by the query

        Returns
        -------
        `numpy.recarray` with OpSim records. The column names may be obtained as
        res.dtype.names

        .. notes:: The `limit` argument should only be used if a small example
        is required. The angle ranges in the argument should be specified in degrees.
        """
        query = self.baseQuery + ' FROM SUMMARY'

        nConstraints = 0 # the number of constraints in this query

        for column in self.columnMapping:
            value = eval(column[0])
            if value is not None:
                if nConstraints > 0:
                    query += ' AND'
                else:
                    query += ' WHERE '

                if isinstance(value, tuple):
                    if len(value)>2:
                        raise RuntimeError('Cannot pass a tuple longer than 2 elements '+
                                           'to getObservationMetaData: %s is len %d'
                                           % (column[0], len(value)))

                    # perform any necessary coordinate transformations
                    if column[4] is not None:
                        vmin = column[4](value[0])
                        vmax = column[4](value[1])
                    else:
                        vmin = value[0]
                        vmax = value[1]

                    query += ' %s > %s AND %s < %s' % \
                             (column[1], vmin, column[1], vmax)
                else:
                    # perform any necessary coordinate transformations
                    if column[4] is not None:
                        vv = column[4](value)
                    else:
                        vv = value
                    query += ' %s == %s' % (column[1], vv)

                nConstraints += 1

        query += ' GROUP BY expMJD'

        if limit is not None:
            query += ' LIMIT %d' % limit

        if nConstraints == 0 and limit is None:
            raise RuntimeError('You did not specify any contraints on your query;' +
                               ' you will just return ObservationMetaData for all poitnings')

        results = self.opsimdb.execute_arbitrary(query, dtype=self.dtype)
        return results


    @staticmethod
    def ObservationMetaDataFromPointing(OpSimPointingRecord, columnMap,
                                       OpSimColumns=None,
                                       boundLength=1.75, boundType='circle'):
        """
        Return instance of ObservationMetaData for an OpSim Pointing record
        from OpSim.

        Parameters
        ----------
        OpSimPointingRecord : Dictionary, mandatory
            Dictionary of values with keys corresponding to certain columns of
            the Summary table in the OpSim database. The minimal list of keys
            required for catsim to work is 'fiveSigmaDepth',
            'filtSkyBrightness', and at least one of ('finSeeing', 'FWHMeff').
            More keys defined in columnMap may be necessary for PhoSim to work.
        columnMap : A mapping of columns with names exepected by phosim for
            OpSim columns. This can be obtained from
            ObservationMetaDataGenerator.OpSimColumnMap()
        OpSimColumns : tuple of strings, optional, defaults to None
            The columns corresponding to the OpSim records. If None, attempts
            to obtain these from the OpSimRecord as OpSimRecord.dtype.names
        boundType : {'circle', 'box'}, optional, defaults to 'circle'
            Shape of the observation
        boundLength : scalar float, optional, defaults to 1.75
            'characteristic size' of observation field, in units of degrees.
            For boundType='circle', this is a radius, for boundType='box', this
            is a size of the box
        """

        pointing = OpSimPointingRecord
        pointing_column_names = pointing.dtype.names
        # Decide what is the name of the column in the OpSim database
        # corresponding to the Seeing. For older OpSim outputs, this is
        # 'finSeeing'. For later OpSim outputs this is 'FWHMeff'
        if OpSimColumns is None:
            OpSimColumns = pointing_column_names

        if 'FWHMeff' in OpSimColumns:
            _seeingColumn = 'FWHMeff'
        elif 'finSeeing' in OpSimColumns:
            _seeingColumn = 'finSeeing'
        else:
            raise ValueError('OpSimPointingRecord has to contain seeing values'
                             "either as 'FWHMeff', or as 'finSeeing'\n")
        ActiveColumns = list(column for column in columnMap
                             if column[1] in OpSimColumns)
        # convert list of tuples of the form (Name, (value, dtype)) to
        # an ordered Dict
        phosimDict = dict([(col[2], pointing[col[1]])
                           for col in columnMap
                           if (col[1] in pointing_column_names and
                               col[1] not in
                               ('fieldRA', 'fieldDec', 'expMJD',
                                'filter', 'fiveSigmaDepth',
                                'filtSkyBrightness', _seeingColumn,
                                'rotSkyPos', 'altitude', 'azimuth',
                                'airmass'))])

        for col in columnMap:
            if col[2] in phosimDict and col[5] is not None:
                phosimDict[col[2]] = col[5](phosimDict[col[2]])

        obs = ObservationMetaData(pointingRA=numpy.degrees(pointing['fieldRA']),
                                  pointingDec=numpy.degrees(pointing['fieldDec']),
                                  mjd=pointing['expMJD'],
                                  rotSkyPos=numpy.degrees(pointing['rotSkyPos']),
                                  bandpassName=pointing['filter'],
                                  boundType=boundType,
                                  boundLength=boundLength)

        if 'fiveSigmaDepth' in pointing_column_names:
            obs.m5 = pointing['fiveSigmaDepth']
        if 'filtSkyBrightness' in pointing_column_names:
            obs.skyBrightness = pointing['filtSkyBrightness']
        if _seeingColumn in pointing_column_names:
            obs.seeing = pointing[_seeingColumn]

        obs.phoSimMetaData = phosimDict

        return obs

    @staticmethod 
    def ObservationMetaDataFromPointingArray(OpSimPointingRecords,
                                                OpSimColumns=None,
                                                boundLength=1.75,
                                                boundType='circle'):
        """
        Static method to get a list of instances of ObservationMetaData
        corresponding to the records in `numpy.recarray`, where it uses
        the dtypes of the recArray for ObservationMetaData attributes that
        require the dtype.
    
        Parameters
        ----------
        OpSimPointingRecords : `numpy.recarray` of OpSim Records
        OpSimColumns : a tuple of strings, optional, defaults to None
            tuple of column Names of the data in the `numpy.recarray`. If
            None, these names are extracted from the recarray.
        boundType : {'circle' or 'box'}
            denotes the shape of the pointing
        boundLength : float, optional, defaults to 1.75
            the bound length of the Pointing in units of degrees. For boundType
            'box', this is the length of the side of the square box. For boundType
            'circle' this is the radius.
        """
    
        if OpSimColumns is None:
            OpSimColumns = OpSimPointingRecords.dtype.names
    
        # Find out what the Seeing Variable is called in these OpSim records
        seeingVar = None
        matches = 0 
        for var in ['finSeeing', 'FWHMeff']:
            if var in OpSimColumns:
                seeingVar = var
                matches += 1
        if matches in [0, 2]:
            raise ValueError('finSeeing or FWHMeff not in OpSimColumn\n')
    
        columnMapping = ObservationMetaDataGenerator.OpSimColumnMap(seeingVar)
    
        out = list(ObservationMetaDataGenerator.ObservationMetaDataFromPointing(OpSimPointingRecord,
                                                        columnMap=columnMapping,
                                                        OpSimColumns=OpSimColumns,
                                                        boundLength=boundLength,
                                                        boundType=boundType)
                   for OpSimPointingRecord in OpSimPointingRecords)
    
        return out

    def getObservationMetaData(self, obsHistID=None, expDate=None, fieldRA=None, fieldDec=None,
                               moonRA=None, moonDec=None, rotSkyPos=None, telescopeFilter=None,
                               rawSeeing=None, seeing=None, sunAlt=None, moonAlt=None, dist2Moon=None,
                               moonPhase=None, expMJD=None, altitude=None, azimuth=None,
                               visitExpTime=None, airmass=None, skyBrightness=None,
                               m5=None, boundType='circle', boundLength=0.1, limit=None):

        """
        This method will query the OpSim database summary table according to user-specified
        constraints and return a list of of ObservationMetaData instantiations consistent
        with those constraints.

        @param [in] limit is an integer denoting the maximum number of ObservationMetaData to
        be returned

        @param [in] boundType is the boundType of the ObservationMetaData to be returned
        (see documentation in sims_catalogs_generation/../db/spatialBounds.py for more
        details)

        @param [in] boundLength is the boundLength of the ObservationMetaData to be
        returned (in degrees; see documentation in
        sims_catalogs_generation/../db/spatialBounds.py for more details)

        All other input parameters are constraints to be placed on the SQL query of the
        opsim output db.  These contraints can either be tuples of the form (min, max)
        or an exact value the user wants returned.

        Parameters that can be constrained are:

        @param [in] fieldRA in degrees
        @param [in] fieldDec in degrees
        @param [in] altitude in degrees
        @param [in] azimuth in degrees

        @param [in] moonRA in degrees
        @param [in] moonDec in degrees
        @param [in] moonAlt in degrees
        @param [in] moonPhase (a value from 1 to 100 indicating how much of the moon is illuminated)
        @param [in] dist2Moon the distance between the telescope pointing and the moon in degrees

        @param [in] sunAlt in degrees

        @param [in[ rotSkyPos (the angle of the sky with respect to the camera coordinate system) in degrees
        @param [in] telescopeFilter a string that is one of u,g,r,i,z,y

        @param [in] airmass
        @param [in] rawSeeing (this is an idealized seeing at zenith at 500nm in arcseconds)
        @param [in] seeing (this is the OpSim column 'FWHMeff' or 'finSeeing' [deprecated] in arcseconds)

        @param [in] visitExpTime the exposure time in seconds
        @param [in] obsHistID the integer used by OpSim to label pointings
        @param [in] expDate is the date of the exposure (units????)
        @param [in] expMJD is the MJD of the exposure
        @param [in] m5 is the five sigma depth of the observation
        @param [in] skyBrightness
        """

        OpSimPointingRecords = self.getOpSimRecords(obsHistID=obsHistID,
                                                    expDate=expDate,
                                                    fieldRA=fieldRA,
                                                    fieldDec=fieldDec,
                                                    moonRA=moonRA,
                                                    moonDec=moonDec,
                                                    rotSkyPos=rotSkyPos,
                                                    telescopeFilter=telescopeFilter,
                                                    rawSeeing=rawSeeing,
                                                    seeing=seeing,
                                                    sunAlt=sunAlt,
                                                    moonAlt=moonAlt,
                                                    dist2Moon=dist2Moon,
                                                    moonPhase=moonPhase,
                                                    expMJD=expMJD,
                                                    altitude=altitude,
                                                    azimuth=azimuth,
                                                    visitExpTime=visitExpTime,
                                                    airmass=airmass,
                                                    skyBrightness=skyBrightness,
                                                    m5=m5, boundType=boundType,
                                                    boundLength=boundLength,
                                                    limit=limit)

        output = self.ObservationMetaDataFromPointingArray(OpSimPointingRecords,
                                                           OpSimColumns=None,
                                                           boundType=boundType,
                                                           boundLength=boundLength)
        return output

        # OpSimColumns = OpSimPointingRecords.dtype.names
        # convert the results into ObservationMetaData instantiations
        # out = list(self.ObservationMetaDataFromPointing(OpSimPointingRecord,
        #                                                 columnMap=self.columnMapping,
        #                                                 OpSimColumns=OpSimColumns,
        #                                                 boundLength=boundLength,
        #                                                 boundType=boundType)
        #          for OpSimPointingRecord in OpSimPointingRecords)
        # return out
