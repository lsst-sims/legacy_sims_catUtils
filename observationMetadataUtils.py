import math
import numpy
from collections import OrderedDict

def equationOfEquinoxes(d):
    l = 280.47 + 0.98565*d
    omega = 125.04 - 0.052954*d
    deltaPsi = -0.000319*math.sin(omega) - 0.000024*math.sin(2*l)
    epsilon = 23.4393 - 0.0000004*d
    return  deltaPsi*math.cos(epsilon)

def calcGmstGast(mjd):
    #From http://aa.usno.navy.mil/faq/docs/GAST.php Nov. 9 2013
    mjdConv = 2400000.5
    jd2000 = 2451545.0
    mjd_o = math.floor(mjd)
    jd = mjd + mjdConv
    jd_o = mjd_o + mjdConv
    h = 24.*(jd-jd_o)
    d = jd - jd2000
    d_o = jd_o - jd2000
    t = d/36525.
    gmst = 6.697374558 + 0.06570982441908*d_o + 1.00273790935*h + 0.000026*t**2
    gast = gmst + equationOfEquinoxes(d)
    gmst %= 24.
    gast %= 24.
    return {'GMST':gmst, 'GAST':gast}

def calcLmstLast(mjd, longRad):
    longDeg = math.degrees(longRad)
    longDeg %= 360.
    if longDeg > 180.:
        longDeg -= 360.
    hrs = longDeg/15.
    gmstgast = calcGmstGast(mjd)
    lmst = gmstgast['GMST']+hrs
    last = gmstgast['GAST']+hrs
    lmst %= 24.
    last %= 24.
    return {'LMST':lmst, 'LAST':last}

def raDecToAltAz(raRad, decRad, longRad, latRad, mjd):
    lst = calcLmstLast(mjd, longRad)
    last = lst['LAST']
    haRad = math.radians(last*15.) - raRad
    altRad = math.asin(math.sin(decRad)*math.sin(latRad)+math.cos(decRad)*math.cos(latRad)*math.cos(haRad))
    azRad = math.acos((math.sin(decRad) - math.sin(altRad)*math.sin(latRad))/(math.cos(altRad)*math.cos(latRad)))
    if math.sin(haRad) >= 0:
        azRad = 2.*math.pi-azRad
    return altRad, azRad

def altAzToRaDec(altRad, azRad, longRad, latRad, mjd):
    lst = calcLmstLast(mjd, longRad)
    last = lst['LAST']
    decRad = math.asin(math.sin(latRad)*math.sin(altRad)+ math.cos(latRad)*math.cos(altRad)*math.cos(azRad))
    haRad = math.acos((math.sin(altRad) - math.sin(decRad)*math.sin(latRad))/(math.cos(decRad)*math.cos(latRad)))
    raRad = math.radians(last*15.) - haRad
    return raRad, decRad

def calcPa(azRad, decRad, latRad):
    """
    Calculate the Parallactic angle 
    azRad is the azimuth of the object assuming OpSim conventions (radians)
    latRad is the latitude of the observatory (radians)
    decRad is the declination of the object (radians)
    """
    try:
        paRad = math.asin(math.sin(azRad)*math.cos(latRad)/math.cos(decRad))
    except ValueError, e:
        if not math.fabs(decRad) > math.fabs(latRad):
            raise ValueError("The point is circumpolar but the Azimuth is not valid: Az=%.2f"%(math.degrees(azRad)))
        else:
            raise e
    return paRad

def getRotSkyPos(azRad, decRad, latRad, rotTelRad):
    """
    azRad is the azimuth of the object assuming opSim conventions (radians)
    decRad is the declination of the object (radians)
    latRad is the latitude of the observatory (radians)
    rotTelRad is the angle of the camera rotator assuming OpSim
    conventions (radians)
    """
    paRad = calcPa(azRad, decRad, latRad)
    return (rotTelRad - paRad + math.pi)%(2.*math.pi)

def getRotTelPos(azRad, decRad, latRad, rotSkyRad):
    """
    azRad is the azimuth of the object assuming opSim conventions (radians)
    decRad is the declination of the object (radians)
    latRad is the latitude of the observatory (radians)
    rotSkyRad is the angle of the field of view relative to the South pole given
    a rotator angle in OpSim conventions (radians)
    """
    paRad = calcPa(azRad, decRad, latRad)
    return (rotSkyRad + paRad - math.pi)%(2.*math.pi)

def haversine(long1, lat1, long2, lat2):
    #From http://en.wikipedia.org/wiki/Haversine_formula
    t1 = math.sin(lat2/2.-lat1/2.)**2
    t2 = math.cos(lat1)*math.cos(lat2)*math.sin(long2/2. - long1/2.)**2
    return 2*math.asin(math.sqrt(t1 + t2))

def calcObsDefaults(raRad, decRad, altRad, azRad, rotTelRad, mjd, band, longRad, latRad):
    obsMd = {}
    #Defaults
    moonra, moondec = altAzToRaDec(-math.pi/2., 0., longRad, latRad, mjd)
    sunalt = -math.pi/2.
    moonalt = -math.pi/2.
    dist2moon = haversine(moonra, moondec, raRad, decRad)
    obsMd['Opsim_moonra'] = moonra
    obsMd['Opsim_moondec'] = moondec
    obsMd['Opsim_sunalt'] = sunalt
    obsMd['Opsim_moonalt'] = moonalt
    obsMd['Opsim_dist2moon'] = dist2moon

    rotSkyPos = getRotSkyPos(azRad, decRad, latRad, rotTelRad)
    obsMd['Opsim_filter'] = band
    obsMd['Unrefracted_RA'] = raRad
    obsMd['Unrefracted_Dec'] = decRad
    obsMd['Opsim_rotskypos'] = rotSkyPos
    obsMd['Opsim_rottelpos'] = rotTelRad
    obsMd['Opsim_altitude'] = altRad
    obsMd['Opsim_azimuth'] = azRad
    return obsMd

def makeObservationMetadata(metaData):
    return OrderedDict([(k,(metaData[k], numpy.asarray(metaData[k]).dtype))
                         for k in metaData])

def makeObsParamsAzAltTel(azRad, altRad, mjd, band, rotTelRad=0., longRad=-1.2320792, latRad=-0.517781017, **kwargs):
    '''
    Calculate a minimal set of observing parameters give the ra, dec, and time of the observation.
    altRad -- Altitude of the boresite of the observation in radians
    azRad -- Azimuth of the boresite of the observation in radians
    mjd -- MJD of the observation
    band -- bandpass of the observation e.g. 'r'
    rotTelRad -- Rotation of the camera relative to the telescope in radians Default=0.
    longRad -- Longitude of the observatory in radians Default=-1.2320792
    latRad -- Latitude of the observatory in radians Default=-0.517781017
    **kwargs -- The kwargs will be put in the returned dictionary overriding the default value if it exists
    '''

    raRad, decRad = altAzToRaDec(altRad, azRad, longRad, latRad, mjd)
    obsMd = calcObsDefaults(raRad, decRad, altRad, azRad, rotTelRad, mjd, band, longRad, latRad)
    obsMd.update(kwargs)
    return makeObservationMetadata(obsMd)

def makeObsParamsAzAltSky(azRad, altRad, mjd, band, rotSkyRad=math.pi, longRad=-1.2320792, latRad=-0.517781017, **kwargs):
    '''
    Calculate a minimal set of observing parameters give the ra, dec, and time of the observation.
    altRad -- Altitude of the boresite of the observation in radians
    azRad -- Azimuth of the boresite of the observation in radians
    mjd -- MJD of the observation
    band -- bandpass of the observation e.g. 'r'
    rotTelRad -- Rotation of the field of view relative to the North pole in radians Default=0.
    longRad -- Longitude of the observatory in radians Default=-1.2320792
    latRad -- Latitude of the observatory in radians Default=-0.517781017
    **kwargs -- The kwargs will be put in the returned dictionary overriding the default value if it exists
    '''
    raRad, decRad = altAzToRaDec(altRad, azRad, longRad, latRad, mjd)
    rotTelRad = getRotTelPos(azRad, decRad, latRad, rotSkyRad)
    return makeObsParamsAzAltTel(azRad, altRad, mjd, band, rotTelRad=rotTelRad, longRad=longRad, latRad=latRad, **kwargs)


def makeObsParamsRaDecTel(raRad, decRad, mjd, band, rotTelRad=0., longRad=-1.2320792, latRad=-0.517781017, **kwargs):
    '''
    Calculate a minimal set of observing parameters give the ra, dec, and time of the observation.
    raRad -- RA of the boresite of the observation in radians
    decRad -- Dec of the boresite of the observation in radians
    mjd -- MJD of the observation
    band -- bandpass of the observation e.g. 'r'
    rotTelRad -- Rotation of the camera relative to the telescope in radians Default=0.
    longRad -- Longitude of the observatory in radians Default=-1.2320792
    latRad -- Latitude of the observatory in radians Default=-0.517781017
    **kwargs -- The kwargs will be put in the returned dictionary overriding the default value if it exists
    '''
    altRad, azRad = raDecToAltAz(raRad, decRad, longRad, latRad, mjd)
    obsMd = calcObsDefaults(raRad, decRad, altRad, azRad, rotTelRad, mjd, band, longRad, latRad)
    obsMd.update(kwargs)
    return makeObservationMetadata(obsMd)

def makeObsParamsRaDecSky(raRad, decRad, mjd, band, rotSkyRad=math.pi, longRad=-1.2320792, latRad=-0.517781017, **kwargs):
    '''
    Calculate a minimal set of observing parameters give the ra, dec, and time of the observation.
    raRad -- RA of the boresite of the observation in radians
    decRad -- Dec of the boresite of the observation in radians
    mjd -- MJD of the observation
    band -- bandpass of the observation e.g. 'r'
    rotSkyRad -- Rotation of the field of view relative to the North pole in radians Default=0.
    longRad -- Longitude of the observatory in radians Default=-1.2320792
    latRad -- Latitude of the observatory in radians Default=-0.517781017
    **kwargs -- The kwargs will be put in the returned dictionary overriding the default value if it exists
    '''
    altRad, azRad = raDecToAltAz(raRad, decRad, longRad, latRad, mjd)
    rotTelRad = getRotTelPos(azRad, decRad, latRad, rotSkyRad)
    return makeObsParamsRaDecTel(raRad, decRad, mjd, band, rotTelRad=rotTelRad, longRad=longRad, latRad=latRad, **kwargs)
