import warnings
from sqlalchemy.orm import scoped_session, sessionmaker, mapper
from sqlalchemy.sql import expression
from sqlalchemy import create_engine
from sqlalchemy import ThreadLocalMetaData
import sqlalchemy.databases as sd 
from sqlalchemy import func
from sqlalchemy import schema
from sqlalchemy import MetaData
from sqlalchemy import Table
from sqlalchemy.ext.sqlsoup import SqlSoup
from sqlalchemy import exc as sa_exc

sessions = {}
warnings.simplefilter("ignore", category=sa_exc.SAWarning)
a_engine =create_engine("mssql+pymssql://LSST-2:L$$TUser@fatboy.npl.washington.edu:1433/LSST",
        echo=False)
a_session = scoped_session(sessionmaker(autoflush=True, 
    bind=a_engine))
a_metadata = MetaData()
a_metadata.bind = a_engine
opsim = Table('output_opsim3_61', a_metadata, autoload=True)
allstar = Table('starsALL_forceseek', a_metadata, autoload=True)
star = Table('starsMSRGB_forceseek', a_metadata, autoload=True)
starcirc = Table('starsMSRGBp80m10_forceseek', a_metadata, autoload=True)
wd = Table('starsWD_forceseek', a_metadata, autoload=True)
rrly = Table('starsRRLy_forceseek', a_metadata, autoload=True)
bhb = Table('starsBHB_forceseek', a_metadata, autoload=True)
lens = Table('lens', a_metadata, autoload=True)
image = Table('image', a_metadata, autoload=True)
eb = Table('ebstars', a_metadata, autoload=True)
cepheid = Table('cepheidstars', a_metadata, autoload=True)
astromeggs = Table('AstromEasterEggs', a_metadata, autoload=True)
dwarfgal = Table('dwarfGalaxies', a_metadata, autoload=True)
db = SqlSoup(a_metadata)
Star = db.map(star, primary_key=[star.c.simobjid])
StarCirc = db.map(starcirc, primary_key=[starcirc.c.simobjid])
AllStar = db.map(allstar, primary_key=[allstar.c.simobjid])
AstromEggs = db.map(astromeggs)
OpSim3_61 = db.map(opsim)
Wd = db.map(wd, primary_key=[wd.c.simobjid])
BHB = db.map(bhb, primary_key=[bhb.c.simobjid])
RRLy = db.map(rrly, primary_key=[rrly.c.simobjid])
LENS = db.map(lens, primary_key=[lens.c.id])
IMAGE = db.map(image, primary_key=[image.c.id])
EBSTARS = db.map(eb, primary_key=[eb.c.simobjid])
CEPHEIDSTARS = db.map(cepheid, primary_key=[cepheid.c.simobjid])
DWARFGALAXYSTARS = db.map(dwarfgal)

sessions['MSSQL'] = a_session
def initGalaxy(ra, dec, radiusdeg, columns, constraint=None):
    if constraint is not None:
        query = a_session.execute("EXECUTE [LSST].[dbo].[GalaxySearchSpecColsConstraint]\
            @RaSearch = %f, @DecSearch = %f, @apertureRadius = %f,\
            @ColumnNames = '%s', @WhereClause =\
            '%s'"%(ra,dec,radiusdeg*60.,columns,constraint))
    else:
        query = a_session.execute("EXECUTE [LSST].[dbo].[GalaxySearchSpecColsConstraint]\
            @RaSearch = %f, @DecSearch = %f, @apertureRadius = %f,\
            @ColumnNames = '%s'"%(ra,dec,radiusdeg*60.,columns))
    coldesc = []
    for k in query.keys():
        coldesc.append({"name":k})
    return query, coldesc

def initSSM(ra, dec, radiusdeg, expmjd, columns, constraint=None):
    if constraint is not None:
        query = a_session.execute("select %s from [LSST].[dbo].[fSSMAll](\
            %f, %f, %f,%f) where %s"%(columns, expmjd, ra, dec, radiusdeg*60., constraint))
    else:
       query = a_session.execute("select %s from [LSST].[dbo].[fSSMAll](\
            %f, %f, %f,%f)"%(columns, expmjd, ra, dec, radiusdeg*60.))
    coldesc = []
    for k in query.keys():
        coldesc.append({"name":k})
    return query, coldesc
