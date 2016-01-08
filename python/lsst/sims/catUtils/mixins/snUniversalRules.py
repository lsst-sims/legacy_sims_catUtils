import numpy as np

__all__ = ['SNUniverse']

class SNUniverse(object):
    """
    Mixin Class for `lsst.sims.catalogs.measures.instances.InstanceCatalog` 
    providing methods for distributing SN model parameters
    conditioned  on the parameters of the hostgalaxy from catsim Catalogs. 
    In particular, this is meant to be used to extend another mixin
    `lsst.sims.catUtils.mixins.sncat.py`

    The base class must have the following attributes:

    numobjs : integer, number of objects in the catalog. 
    badvalues : a value to be reported when calculations fail, typically
        `np.nan`
    suppressHighzSN : Bool, if true drop any SN with redshift higher than
        a threshold `self.maxz`
    maxz : float, threshold redshift value, so that SN with z > self.maxz are
        dropped if self.suppressHighzSN is True
    midSurveyTime : float, MJD at the middle of the survey
    mjdobs : float, MJD of the date of observation of the instance catalog
    maxTimeSNVisible : float, units of days, Time beyond which SN is assumed
        to be invisible. (see definition of suppressDimSN)
    suppressDimSN : bool, if true drop observations of SN with peak at t0, if
        abs(mjdobs - t0) > self.maxTimeVisible

    """

    @property
    def snFrequency(self):
        """
        Frequency with which SN occur in the host galaxy in units of days^{-1}
        """
        if not hasattr(self, '_snFrequency'):
            hundredyear = 100. * 365.
            self._snFrequency = 1.0 / hundredyear
        return self._snFrequency

    @snFrequency.setter
    def snFrequency(self, value):
        """
        set the value of snFrequency

        Parameters
        ----------
        value : float, mandatory
            value to which to set snFrequency to (units of 1/days)
        """
        self._snFrequency = value
        return self._snFrequency
            
    @property
    def midSurveyTime(self):
        '''
        The time at the middle of the survey: ie. at the 5 year period.


        .. note: Changing this should not change the statistical
        properties of the survey, but will change the exact SN we find.
        '''
        if not hasattr(self, '_midSurveyTime'):
            midSurveyTime_default = 570000.0
            self._midSurveyTime = midSurveyTime_default
        return self._midSurveyTime

    @midSurveyTime.setter
    def midSurveyTime(self, mymidSurveyTime):
        '''
        set the value of `self.midSurveyTime` of the catalog

        Parameters
        ----------
        mymidSurveyTime : float, mandatory
            Value to set the midSurveyTime attribute to
        '''
        # if suppressDimSN is None:
        #    self._suppressDimSN = True
        # else:
        self._midSurveyTime = mymidSurveyTime
        return self._midSurveyTime



    def SNCoordinatesFromHost(self, hostra, hostdec, hostz):
        '''
        Distribution of SN coordinates and velocities given a set of host
        coordinates and velocities.
        '''
        numhosts = self.numobjs 

        sndec = hostdec
        snra = hostra
        snz = hostz

        snvra = np.zeros(numhosts)
        snvdec = np.zeros(numhosts)
        snvr = np.zeros(numhosts)

        if self.suppressHighzSN:
            snz = np.where(snz > self.maxz, np.nan, snz)
    
        return snra, sndec, snz , snvra, snvdec, snvr 

    def SNparamDistfromHost(self, hostz, hostid, hostmu):
        '''
        Distribution of SN model parameters given their hosts
        '''
        vals = np.zeros(shape=(self.numobjs, 4))

        for i, v in enumerate(vals):
            vals[i, :] = self.drawSNParams(hostid[i], hostmu[i])
        
        return vals
    def drawSNParams(self, hostid, hostmu):
        """
        return the SALT2 parameters c, x1, x0, t0 for a SN with a particular
        hostid/seed, and distance modulus value in Mpc

        Parameters
        ----------
        hostid: int, mandatory

        hostmu: float, mandatory
        """
        hostid = hostid % 4294967295
        np.random.seed(hostid)
        t0val = self.drawFromT0Dist() 
        if t0val is self.badvalues:
            return [self.badvalues]*4
        else:
            cval = self.drawFromcDist()
            x1val = self.drawFromx1Dist()
            x0val = self.drawFromX0Dist(x1val, cval, hostmu=hostmu)
        return [cval, x1val, x0val, t0val]



                    
            
    def drawFromx1Dist(self, **hostParams):
        """
        """
        return np.random.normal(0. , 1.0) 

    def drawFromcDist(self, **hostParams):
        """
        """
        return np.random.normal(0. , 0.1) 
    
    def drawFromX0Dist(self, x1val, cval, hostmu, **hostParams):
        """
        """
        import snObject

        # First draw an absolute BessellB magnitude for SN
        mabs = np.random.normal(-19.3, 0.3)
        mag = mabs + hostmu

        sn = snObject.SNObject()
        sn.set(x1=x1val, c=cval)
        sn.source.set_peakmag(mag, band='bessellb', magsys='ab')
        x0val = sn.get('x0')
        
        return x0val

    def drawFromT0Dist(self, **hostParams):
        '''
        Distribution function of the time of peak of SN

        '''

        # Will not use hostid for now, but this is there so that 
        # later on one could obtain z, hostmass etc. This is useful to obtain
        # z and host dependent SN Frequency
        hundredyear = 1.0 / self.snFrequency
        t0val = np.random.uniform(-hundredyear / 2.0 + self.midSurveyTime, 
                           hundredyear / 2.0 + self.midSurveyTime)

        if self.suppressDimSN:

            # SN far away from peak will have zero flux, but writing zero flux
            # objects to catalogs will mean each instance catalog will have
            # too many objects. So, this variable will always stay True (except
            # for testing purposes), while the variable that will be changed
            # is maxTimeSNVisible

            if np.abs(t0val - self.mjdobs) > self.maxTimeSNVisible:
      return t0val

