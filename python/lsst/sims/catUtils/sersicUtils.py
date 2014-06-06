from scipy.special import gammainc
from scipy.interpolate import interp1d
import numpy

class SersicUtils(object):
    
    def __init__(self, comp_n, min_log_s=-6, max_log_s=6, dlog_s=0.001):
        self.minprob = []
        self.maxprob = []
        self.interparr = []
        log_s = numpy.arange(min_log_s, max_log_s, dlog_s)
        tinterp, minp, maxp = self._getInterpProbs(comp_n, log_s)
        self.pinterp = tinterp
        self.minprob = minp
        self.maxprob = maxp

    def _makeArray(self, arr):
        if hasattr(arr, "__iter__"):
            return numpy.array(arr)
        else:
            return numpy.array([arr,])

    def _getInterpProbs(self, n, log_s):
        #These values come from:
        #http://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
        b_n = 1.9992*n-0.3271
        x = b_n*numpy.power(10,log_s)**(1./n)
        probs = gammainc(2*n, x)
        return interp1d(probs, log_s), numpy.min(probs), numpy.max(probs)

    def getR(self, Re):
        Re = self._makeArray(Re)
        rands = numpy.random.random(len(Re))
        idxs = numpy.where(rands < self.minprob)[0]
        rands[idxs] = self.minprob
        idxs = numpy.where(rands > self.maxprob)[0]
        rands[idxs] = self.maxprob
        svals = numpy.power(10., self.pinterp(rands))
        return Re*svals

    def getRe(self, a, b, angles):
        a = self._makeArray(a)
        b = self._makeArray(b)
        angles = self._makeArray(angles)
        return (a*b)/numpy.sqrt((b*numpy.cos(angles))**2 + (a*numpy.sin(angles))**2)

    def getXY(self, R, angles):
        x = R*numpy.cos(angles)
        y = R*numpy.sin(angles)
        return x, y


    def getAngles(self, a, b, num):
        #from wolfram:
        #http://mathworld.wolfram.com/Ellipse.html
        e = b/a
        tvals = numpy.random.random(num)*2.*numpy.pi
        theta = numpy.arctan(e*numpy.tan(tvals))%numpy.pi
        idx = numpy.where(tvals > numpy.pi)[0]
        theta[idx] += numpy.pi
        return theta

    def calcDistribution(self, a, b, pa, num):
        angles = self.getAngles(a, b, num)
        Res = self.getRe(a, b, angles)
        Rs = self.getR(Res)
        return self.getXY(Rs, angles+pa)
    
    @staticmethod
    def getEllip(x, y):
        Ixx = numpy.sum(x*x)/len(x)
        Iyy = numpy.sum(y*y)/len(y)
        Ixy = numpy.sum(x*y)/len(x)
        rs = numpy.sqrt(x*x + y*y)
        rs = numpy.sort(rs)
        lrs = len(rs)
        if lrs%2 == 0:
            r_hl = (rs[lrs/2] + rs[lrs/2-1])/2.
        else:
            r_hl = rs[lrs/2]
        #As defined in Williams et al. 1996
        r1 = numpy.sum(numpy.sqrt(x*x + y*y))/len(x)
        e1 = (Ixx - Iyy)/(Ixx + Iyy + 2*numpy.sqrt(Ixx*Iyy - Ixy*Ixy))
        e2 = 2*Ixy/(Ixx + Iyy + 2*numpy.sqrt(Ixx*Iyy - Ixy*Ixy))
        return e1, e2, Ixx, Iyy, Ixy, r1, r_hl
    
    @staticmethod
    def getNphots(ntot, magBulgeNorm, magDiskNorm):
         if numpy.isnan(magBulgeNorm) or magBulgeNorm < 5.:
             return 0, ntot
         elif numpy.isnan(magDiskNorm) or magDiskNorm < 5.:
             return ntot, 0
         else:
            a = numpy.power(10.,(magBulgeNorm - magDiskNorm)/-2.5)
            if numpy.isinf(a):
                return 0, ntot
            b = a/(1.+a)
            return ntot*b, ntot*(1.-b)
