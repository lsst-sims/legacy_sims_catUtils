import math
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.utils import ObservationMetaData, SpatialBounds
#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catUtils.exampleCatalogDefinitions import RefCatalogGalaxyBase, PhoSimCatalogPoint,\
                                                         PhoSimCatalogZPoint, PhoSimCatalogSersic2D
from lsst.sims.utils import raDecFromAltAz

def exampleReferenceCatalog():
    """
    This method outputs a reference catalog of galaxies (i.e. a catalog of
    galaxies in which the columns are simply the data stored in the database).

    The catalog class is defined in
    python/lsst/sims/catUtils/exampleCatalogDefinitions/refCatalogExamples.py

    The catalog is output to the file test_reference.dat
    """

    obs_metadata = ObservationMetaData(boundType='circle',pointingRA=0.0,pointingDec=0.0,
                       boundLength=0.01, mjd=53850.0)
    dbobj = CatalogDBObject.from_objid('galaxyBase')

    t = dbobj.getCatalog('ref_catalog_galaxy', obs_metadata=obs_metadata)
    filename = 'test_reference.dat'
    t.write_catalog(filename, chunk_size=10)

def examplePhoSimCatalogs():
    """
    This method outputs several phoSim input files consisting of different
    types of objects(stars, galaxy bulges, galaxy disks, and AGNs)

    The files created are
    catalog_test_msstars.dat
    catalog_test_galaxyDisk.dat
    catalog_test_galaxyBulge.dat
    catalog_test_galaxyAgn.dat

    (versions are also created that end in _chunked.dat; these should have
    the same contents)

    """
    obsMD = bcm.OpSim3_61DBObject()
    obs_metadata_list = obsMD.getObservationMetaData((55.0, -20.0), 1.0, fovRadius=0.1, makeCircBounds=True)
    obs_metadata = obs_metadata_list[0]
    objectDict = {}
    objectDict['testStars'] = {'dbobj':CatalogDBObject.from_objid('msstars'),
                               'constraint':None,
                               'filetype':'phoSim_catalog_POINT',
                               'obsMetadata':obs_metadata}
    objectDict['testGalaxyBulge'] = {'dbobj':CatalogDBObject.from_objid('galaxyBulge'),
                               'constraint':"mass_bulge > 1. and sedname_bulge is not NULL",
                               'filetype':'phoSim_catalog_SERSIC2D',
                               'obsMetadata':obs_metadata}
    objectDict['testGalaxyDisk'] = {'dbobj':CatalogDBObject.from_objid('galaxyDisk'),
                               'constraint':"DiskLSSTg < 20. and sedname_disk is not NULL",
                               'filetype':'phoSim_catalog_SERSIC2D',
                               'obsMetadata':obs_metadata}
    objectDict['testGalaxyAgn'] = {'dbobj':CatalogDBObject.from_objid('galaxyAgn'),
                               'constraint':"sedname_agn is not NULL",
                               'filetype':'phoSim_catalog_ZPOINT',
                               'obsMetadata':obs_metadata}

    for objKey in objectDict.keys():
        dbobj = objectDict[objKey]['dbobj']
        t = dbobj.getCatalog(objectDict[objKey]['filetype'],
                             obs_metadata=objectDict[objKey]['obsMetadata'],
                             constraint=objectDict[objKey]['constraint'])

        t.phoSimHeaderMap = {}  # technically only needed for PhoSim InstanceCatalog classes

        print
        print "These are the required columns from the database:"
        print t.db_required_columns()
        print
        print "These are the columns that will be output to the file:"
        print t.column_outputs
        print

        filename = 'catalog_test_%s.dat'%(dbobj.objid)
        print "querying and writing catalog to %s:" % filename
        t.write_catalog(filename)
        filename = 'catalog_test_%s_chunked.dat'%(dbobj.objid)
        t.write_catalog(filename, chunk_size=10)
        print " - finished"

def examplePhoSimNoOpSim():
    """
    This method outputs phoSim input files based on arbitrary input coordinates
    (rather than an OpSim pointing).

    catalog_test_stars_rd.dat is a file created from a specified RA, Dec pointing

    catalog_test_stars_aa.dat is a file created from a specified Alt, Az pointing
    """
    raDeg= 15.
    decDeg = -30.
    mjd = 51999.75

    obs_metadata_rd = ObservationMetaData(boundType='circle',
                                          boundLength=0.1,
                                          mjd=mjd,
                                          pointingRA=raDeg,
                                          pointingDec=decDeg,
                                          rotSkyPos=22.0,
                                          bandpassName='g')
    az = 220.0
    alt = 79.0
    mjd = 55958.0
    obs_dummy = ObservationMetaData(mjd=mjd)
    ra, dec = raDecFromAltAz(alt, az, obs_dummy)
    obs_metadata_aa = ObservationMetaData(boundType='circle',
                                          boundLength=0.1,
                                          mjd=mjd,
                                          rotSkyPos=22.0,
                                          pointingRA=ra,
                                          pointingDec=dec,
                                          bandpassName='g')

    dbobj = CatalogDBObject.from_objid('msstars')
    t = dbobj.getCatalog('phoSim_catalog_POINT', obs_metadata= obs_metadata_rd)
    t.phoSimHeaderMap = {}
    t.write_catalog('catalog_test_stars_rd.dat')
    t = dbobj.getCatalog('phoSim_catalog_POINT', obs_metadata= obs_metadata_aa)
    t.phoSimHeaderMap = {}
    t.write_catalog('catalog_test_stars_aa.dat')

def exampleAirmass(airmass,ra = 0.0, dec = 0.0, tol = 30.0, radiusDeg = 0.1,
            makeBoxBounds=False, makeCircBounds=True):
    """
    This method will output a catalog of stars based on an OpSim pointing with
    a specific airmass.  It searches OpSim for pointings with the specified airmass
    and RA, Dec a circle on the sky of radius 'tol' centred on 'ra', 'dec'.
    It creates observation metadata out of the first pointing found and
    uses that to construct the catalog.

    The catalog is output to stars_airmass_test.dat
    """

    obsMD=bcm.OpSim3_61DBObject()

    #The code below will query the OpSim data base object created above.
    #The query will be based on a box in RA, Dec and a specific airmass value

    airmassConstraint = "airmass="+str(airmass) #an SQL constraint that the airmass must be equal to
                                                #the passed value

    #convert q into observation meta data for use in a catalog
    obsMetaData_list = obsMD.getObservationMetaData((ra, dec), tol, constraint=airmassConstraint,
                                                    fovRadius=radiusDeg,makeBoxBounds=makeBoxBounds,
                                                    makeCircBounds=makeCircBounds)

    obsMetaData = obsMetaData_list[0]

    #create and output a reference catalog of stars based on our query to opSim
    dbobj = CatalogDBObject.from_objid('allstars')
    catalog = dbobj.getCatalog('ref_catalog_star', obs_metadata = obsMetaData)
    catalog.write_catalog('stars_airmass_test.dat')


if __name__ == '__main__':
    exampleReferenceCatalog()
    examplePhoSimCatalogs()
    examplePhoSimNoOpSim()
    exampleAirmass(1.1)
