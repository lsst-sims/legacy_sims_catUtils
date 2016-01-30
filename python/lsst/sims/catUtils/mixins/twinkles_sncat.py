"""
Redfine some functions in SNIaCatalog for Twinkles so that it can use TwinkSN
"""
from lsst.sims.catalogs.measures.instance import compound
from .sncat import SNIaCatalog

__all__ = ['TwinklesSNCat']

class TwinklesSNCat(SNIaCatalog):

    # The t0 value stored in the database is in terms of MJD - survey start
    # date. survey start date must be stored so that times correspond to
    # ObservationMetaData times
    # This can be reset by the user
    surveyStartDate = 59580. # For Kraken_1042

#    def get_snid(self):
#        """
#        Obtain the snid column stored in the database
#        """
#        # return self.column_by_name(self.refIdCol)
#        return self.column_by_name('snid')

#     @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
#     def get_angularCoordinates(self):
#         """
#         Obtain the positional parameters of the SN stored in the database.
#         """
# 
#         snra = self.column_by_name('snra')
#         sndec = self.column_by_name('sndec')
#         snz = self.column_by_name('redshift')
#         snvra = np.zeros(self.numobjs)
#         snvdec = np.zeros(self.numobjs)
#         snvr = np.zeros(self.numobjs)
#         return ([snra, sndec, snz, snvra, snvdec, snvr])

    @compound('c', 'x1', 'x0', 't0')
    def get_snparams(self):
        """
        Obtain the SN model parameters stored in the database
        """
        snc = self.column_by_name('c')
        snx1 = self.column_by_name('x1')
        snx0 = self.column_by_name('x0')

        # Shift t0 by surveyStartDate
        snt0 = self.column_by_name('t0')
        snt0 += self.surveyStartDate

        return ([snc, snx1, snx0, snt0])
