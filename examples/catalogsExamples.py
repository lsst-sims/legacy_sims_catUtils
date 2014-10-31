import math
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
#The following is to get the object ids in the registry
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.sims.catUtils.exampleCatalogDefinitions import RefCatalogGalaxyBase, PhoSimCatalogPoint,\
                                                         PhoSimCatalogZPoint, PhoSimCatalogSersic2D
from lsst.sims.catalogs.generation.db import makeObsParamsAzAltTel, makeObsParamsRaDecTel, SpatialBounds

def exampleReferenceCatalog():
    """
    This method outputs a reference catalog of galaxies (i.e. a catalog of
    galaxies in which the columns are simply the data stored in the database).

    The catalog class is defined in
    python/lsst/sims/catUtils/exampleCatalogDefinitions/refCatalogExamples.py

    The catalog is output to the file test_reference.dat
    """

    obs_metadata = ObservationMetaData(boundType='circle',unrefractedRA=0.0,unrefractedDec=0.0,
                       boundLength=0.01)
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
    obs_metadata = obsMD.getObservationMetaData(88544919, 0.1, makeCircBounds=True)
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

    #source code for the method below is found in python.lsst.sims.catUtils.observationMetaDataUtils.py
    #
    #basically, it returns a dict of data needed by phoSim to specify a given pointing
    #(ra and dec of the moon, rotation of the sky relative to the telescope, etc.)
    md =  makeObsParamsRaDecTel(math.radians(raDeg), math.radians(decDeg), mjd, 'r')

    obs_metadata_rd = ObservationMetaData(boundType='circle',unrefractedRA=raDeg,unrefractedDec=decDeg,
                                                        boundLength=0.1,
                                                        mjd=mjd,
                                                        bandpassName='r',
                                                        phoSimMetadata=md)
    azRad = math.radians(220.)
    altRad = math.radians(79.)
    md = makeObsParamsAzAltTel(azRad, altRad, mjd, 'r')
    raDeg = math.degrees(md['Unrefracted_RA'][0])
    decDeg = math.degrees(md['Unrefracted_Dec'][0])
    obs_metadata_aa = ObservationMetaData(boundType='circle',unrefractedRA=raDeg,unrefractedDec=decDeg,
                                                        boundLength=0.1,
                                                        mjd=mjd,
                                                        bandpassName='r',
                                                        phoSimMetadata=md)
    dbobj = CatalogDBObject.from_objid('msstars')
    t = dbobj.getCatalog('phoSim_catalog_POINT', obs_metadata= obs_metadata_rd)
    t.write_catalog('catalog_test_stars_rd.dat')
    t = dbobj.getCatalog('phoSim_catalog_POINT', obs_metadata= obs_metadata_aa)
    t.write_catalog('catalog_test_stars_aa.dat')

def exampleAirmass(airmass,ra = 0.0, dec = 0.0, tol = 10.0, radiusDeg = 0.1,
            makeBoxBounds=False, makeCircBounds=True):
    """
    This method will output a catalog of stars based on an OpSim pointing with
    a specific airmass.  It searches OpSim for pointings with the specified airmass
    and RA, Dec within a box bounded by tol (in degrees).  It creates observation meta
    data out of the first pointing found and uses that to construct the catalog.

    The catalog is output to stars_airmass_test.dat
    """

    obsMD=bcm.OpSim3_61DBObject()


    #The code below will query the OpSim data base object created above.
    #The query will be based on a box in RA, Dec and a specific airmass value

    airmassConstraint = "airmass="+str(airmass) #an SQL constraint that the airmass must be equal to
                                                #the passed value
    skyBounds = SpatialBounds.getSpatialBounds('box', ra, dec, tol)
    
    query = obsMD.executeConstrainedQuery(skyBounds, constraint=airmassConstraint)
    
    #convert q into observation meta data for use in a catalog
    obsMetaData = obsMD.getObservationMetaData(query['Opsim_obshistid'][0],radiusDeg,makeBoxBounds=makeBoxBounds,
                   makeCircBounds=makeCircBounds)

    #create and output a reference catalog of stars based on our query to opSim
    dbobj = CatalogDBObject.from_objid('allstars')
    catalog = dbobj.getCatalog('ref_catalog_star', obs_metadata = obsMetaData)
    catalog.write_catalog('stars_airmass_test.dat')


if __name__ == '__main__':
    exampleReferenceCatalog()
    examplePhoSimCatalogs()
    examplePhoSimNoOpSim()
    exampleAirmass(1.1)
