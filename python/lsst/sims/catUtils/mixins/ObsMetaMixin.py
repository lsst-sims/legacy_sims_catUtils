import numpy as np

__all__ = ["ObsMetadataBase"]

class ObsMetadataBase(object):
    """
    This class provides the InstanceCatalog with a way to write out parameters from its ObsMetadata object.
    """
    def _get_nobj(self):
        return len(self.get_objId())

    def get_expMJD(self):
        n_records = self._get_nobj()
        mjd = [self.obs_metadata.mjd] * n_records
        return np.array(mjd)

    def get_rotSkyPos(self):
        n_records = self._get_nobj()
        rotSkyPos = [self.obs_metadata.rotSkyPos] * n_records
        return np.array(rotSkyPos)

    def get_bandpass(self):
        n_records = self._get_nobj()
        bandpass = [self.obs_metadata.bandpass] * n_records
        return np.array(bandpass)

    def get_m5(self):
        n_records = self._get_nobj()
        m5 = [self.obs_metadata.m5] * n_records
        return np.array(m5)

    def get_seeing(self):
        n_records = self._get_nobj()
        seeing = [self.obs_metadata.seeing[self.obs_metadata.bandpass]] * n_records
        return np.array(seeing)

    def get_fieldRA(self):
        n_records = self._get_nobj()
        ra = [self.obs_metadata.unrefractedRA] * n_records
        return np.array(ra)

    def get_fieldDec(self):
        n_records = self._get_nobj()
        dec = [self.obs_metadata.unrefractedDec] * n_records
        return np.array(dec)

    def get_visitExpTime(self):
        n_records = self._get_nobj()
        visitExpTime = [self.obs_metadata.phoSimMetaData['exptime'][0]] * n_records
        return np.array(visitExpTime)
