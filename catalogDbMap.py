from sqlalchemy import schema
import copy
class catalogDbMap (object):
    def __init__(self):
        self.objectTypes = {}
        self.objectTypes['POINT'] = {}
        self.objectTypes['OPSIM361'] = {}
        self.objectTypes['MOVINGPOINT'] = {}
        self.objectTypes['SERSIC2D_BULGE'] = {}
        self.objectTypes['SERSIC2D_DISK'] = {}
        self.objectTypes['SERSIC2D_AGN'] = {}
        self.objectTypes['SERSIC2D'] = {}
        self.objectTypes['STUB'] = {}
        self.objectTypes['POINT']['id'] = 'id'
        self.objectTypes['POINT']['ra'] = 'ra*PI()/180.'
        self.objectTypes['POINT']['dec'] = 'decl*PI()/180.'
        self.objectTypes['POINT']['glon'] = 'gal_l*PI()/180.'
        self.objectTypes['POINT']['glat'] = 'gal_b*PI()/180.'
        self.objectTypes['POINT']['magNorm'] = 'toMag(flux_scale*flux_ratio_from_lc(isvar,t0,%f,timescale,varfluxpeak,\'%c\'))'
        self.objectTypes['POINT']['sedFilename'] = 'sedfilename'
        self.objectTypes['POINT']['redshift'] = '0.0'
        self.objectTypes['POINT']['shearXX'] = '0.0'
        self.objectTypes['POINT']['shearYY'] = '0.0'
        self.objectTypes['POINT']['magnification'] = '0.0'
        self.objectTypes['POINT']['properMotionRa'] = '(mudecl/(1000.*3600.))*PI()/180.'
        self.objectTypes['POINT']['properMotionDec'] = '(mura/(1000.*3600.))*PI()/180.'
        self.objectTypes['POINT']['spatialmodel'] = "'None'"
        self.objectTypes['POINT']['galacticExtinctionModel'] = "'CCM'"
        self.objectTypes['POINT']['galacticAv'] = 'ebv*3.1'
        self.objectTypes['POINT']['galacticRv'] = '3.1::real'
        self.objectTypes['POINT']['internalExtinctionModel'] = "'None'"
        self.objectTypes['POINT']['parallax'] = 'parallax'
        self.objectTypes['POINT']['radialVelocity'] = 'vr'
        self.objectTypes['OPSIM361']['Unrefracted_RA'] = 'fieldradeg'
        self.objectTypes['OPSIM361']['Unrefracted_Dec'] = 'fielddecdeg'
        self.objectTypes['OPSIM361']['Opsim_moonra'] = 'moonra'
        self.objectTypes['OPSIM361']['Opsim_moondec'] = 'moondec'
        self.objectTypes['OPSIM361']['Opsim_rotskypos'] = 'rotskypos'
        self.objectTypes['OPSIM361']['Opsim_rottelpos'] = 'rottelpos'
        self.objectTypes['OPSIM361']['Opsim_filter'] = 'filter'
        self.objectTypes['OPSIM361']['Opsim_rawseeing'] = 'rawseeing'
        self.objectTypes['OPSIM361']['Opsim_sunalt'] = 'sunalt'
        self.objectTypes['OPSIM361']['Opsim_moonalt'] = 'moonalt'
        self.objectTypes['OPSIM361']['Opsim_dist2moon'] = 'dist2moon'
        self.objectTypes['OPSIM361']['Opsim_moonphase'] = 'moonphase'
        self.objectTypes['OPSIM361']['Opsim_obshistid'] = 'obshistid'
        self.objectTypes['OPSIM361']['Opsim_expmjd'] = 'expmjd'
        self.objectTypes['OPSIM361']['SIM_SEED'] = 'expdate'
        self.objectTypes['SERSIC2D']['id'] = 'id'
        self.objectTypes['SERSIC2D']['ra'] = 'ra'
        self.objectTypes['SERSIC2D']['dec'] = 'dec'
        self.objectTypes['SERSIC2D']['glon'] = '0'
        self.objectTypes['SERSIC2D']['glat'] = '0'
        self.objectTypes['SERSIC2D']['magNorm'] = '0.0'
        self.objectTypes['SERSIC2D']['sedFilename'] = "'None'"
        self.objectTypes['SERSIC2D']['redshift'] = 'redshift'
        self.objectTypes['SERSIC2D']['shearXX'] = '0.0'
        self.objectTypes['SERSIC2D']['shearYY'] = '0.0'
        self.objectTypes['SERSIC2D']['magnification'] = '0.0'
        self.objectTypes['SERSIC2D']['properMotionRa'] = '0.0'
        self.objectTypes['SERSIC2D']['properMotionDec'] = '0.0'
        self.objectTypes['SERSIC2D']['spatialmodel'] = "'None'"
        self.objectTypes['SERSIC2D']['majorAxis'] = '0'
        self.objectTypes['SERSIC2D']['minorAxis'] = '0'
        self.objectTypes['SERSIC2D']['positionAngle'] = '0'
        self.objectTypes['SERSIC2D']['galacticExtinctionModel'] = "'CCM'"
        self.objectTypes['SERSIC2D']['galacticAv'] = '0.0'
        self.objectTypes['SERSIC2D']['galacticRv'] = '3.1'
        self.objectTypes['SERSIC2D']['internalExtinctionModel'] = "'CCM'"
        self.objectTypes['SERSIC2D']['internalAv'] = '0.0'
        self.objectTypes['SERSIC2D']['internalRv'] = '3.1'

        self.objectTypes['SERSIC2D_BULGE'] = copy.copy(self.objectTypes['SERSIC2D'])
        self.objectTypes['SERSIC2D_BULGE']['ra'] = 'bra'
        self.objectTypes['SERSIC2D_BULGE']['dec'] = 'bdec'
        self.objectTypes['SERSIC2D_BULGE']['magNorm'] = 'toMag(flux_scale_bulge)'
        self.objectTypes['SERSIC2D_BULGE']['sedFilename'] =\
          'lookupGalaxySedFromId(sedid_bulge)'
        self.objectTypes['SERSIC2D_BULGE']['majorAxis'] = 'a_b'
        self.objectTypes['SERSIC2D_BULGE']['minorAxis'] = 'b_b'
        self.objectTypes['SERSIC2D_BULGE']['positionAngle'] = 'pa_b'
        self.objectTypes['SERSIC2D_BULGE']['internalExtinctionModel'] = "ext_model_b"
        self.objectTypes['SERSIC2D_BULGE']['internalAv'] = 'av_b'
        self.objectTypes['SERSIC2D_BULGE']['internalRv'] = 'r_vb'

        self.objectTypes['SERSIC2D_DISK'] = copy.copy(self.objectTypes['SERSIC2D'])
        self.objectTypes['SERSIC2D_DISK']['ra'] = 'dra'
        self.objectTypes['SERSIC2D_DISK']['dec'] = 'ddec'
        self.objectTypes['SERSIC2D_DISK']['magNorm'] = 'toMag(flux_scale_bulge)'
        self.objectTypes['SERSIC2D_DISK']['sedFilename'] =\
          'lookupGalaxySedFromId(sedid_disk)'
        self.objectTypes['SERSIC2D_DISK']['majorAxis'] = 'a_d'
        self.objectTypes['SERSIC2D_DISK']['minorAxis'] = 'b_d'
        self.objectTypes['SERSIC2D_DISK']['positionAngle'] = 'pa_d'
        self.objectTypes['SERSIC2D_DISK']['internalExtinctionModel'] = "ext_model_d"
        self.objectTypes['SERSIC2D_DISK']['internalAv'] = 'av_d'
        self.objectTypes['SERSIC2D_DISK']['internalRv'] = 'r_vd'

        self.objectTypes['SERSIC2D_AGN'] = copy.copy(self.objectTypes['SERSIC2D'])
        self.objectTypes['SERSIC2D_AGN']['ra'] = 'agnra'
        self.objectTypes['SERSIC2D_AGN']['dec'] = 'agndec'
        self.objectTypes['SERSIC2D_AGN']['magNorm'] = 'toMag(flux_scale_agn)'
        self.objectTypes['SERSIC2D_AGN']['sedFilename'] =\
          'lookupGalaxySedFromId(sedid_agn)'

class physicalObjectMap (object):
    def __init__(self):
        self.objectMap = {}
        self.objectMap['STARS'] = \
            ({'table':'Star','ptype':'POINT','constraint':None},\
             {'table':'Wd','ptype':'POINT','constraint':None})
        self.objectMap['WDSTARS'] = \
            ({'table':'Wd','ptype':'POINT','constraint':None},)
        self.objectMap['GALAXY'] = \
            ({'table':'Galaxy','ptype':'SERSIC2D_BULGE','constraint':'flux_scale_bulge is not NULL'},\
             {'table':'Galaxy','ptype':'SERSIC2D_DISK','constraint':'flux_scale_disk is not NULL'},\
             {'table':'Galaxy','ptype':'SERSIC2D_AGN','constraint':'isagn > 0'})
        self.objectMap['SSM'] = \
            ({'table':'Ephems','ptype':'MOVINGPOINT','constraint':None},)
        self.objectMap['OPSIM361'] = \
            ({'table':'OpSim3_61','ptype':'OPSIM361','constraint':None},)
