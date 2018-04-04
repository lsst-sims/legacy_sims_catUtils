from builtins import object
import numpy as np

__all__ = ["ObsMetadataBase"]

class ObsMetadataBase(object):
    """
    This class provides the InstanceCatalog with a way to write out parameters from its ObsMetadata object.
    """
    def _get_nobj(self):
        return len(self.column_by_name('objId'))

    def get_expMJD(self):
        """
        Returns the expMJD from obs_metadata (assumed to be TAI MJD exposure midpoint).
        """
        n_records = self._get_nobj()
        mjd = [self.obs_metadata.mjd.TAI] * n_records
        return np.array(mjd)

    def get_rotSkyPos(self):
        """
        Returns the rotSkyPos from obs_metadata (in degrees).
        """
        n_records = self._get_nobj()
        rotSkyPos = [self.obs_metadata.rotSkyPos] * n_records
        return np.array(rotSkyPos)

    def get_bandpass(self):
        """
        Returns the bandpass name from obs_metadata.
        """
        n_records = self._get_nobj()
        bandpass = [self.obs_metadata.bandpass] * n_records
        return np.array(bandpass)

    def get_m5(self):
        """
        Returns m5 from obs_metadata.
        """
        n_records = self._get_nobj()
        m5 = [self.obs_metadata.m5[self.obs_metadata.bandpass]] * n_records
        return np.array(m5)

    def get_seeing(self):
        """
        Returns the seeing from obs_metadata (in arcseconds).
        """
        n_records = self._get_nobj()
        seeing = [self.obs_metadata.seeing[self.obs_metadata.bandpass]] * n_records
        return np.array(seeing)

    def get_fieldRA(self):
        """
        Returns the field RA from obs_metadata (in degrees).
        """
        n_records = self._get_nobj()
        ra = [self.obs_metadata.unrefractedRA] * n_records
        return np.array(ra)

    def get_fieldDec(self):
        """
        Returns the field Dec from obs_metadata (in degrees).
        """
        n_records = self._get_nobj()
        dec = [self.obs_metadata.unrefractedDec] * n_records
        return np.array(dec)

    def get_visitExpTime(self):
        """
        Returns the visitExpTime (open shutter time) in seconds.
        """
        n_records = self._get_nobj()
        try:
            visitExpTime = [self.obs_metadata.OpsimMetaData['visitExpTime']] * n_records
        except KeyError:
            try:
                # V4 future proofing.
                visitExpTime = [self.obs_metadata.OpsimMetaData['visitExposureTime']] * n_records
        return np.array(visitExpTime)
