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
#from elixir import *


a_engine =create_engine("mssql://LSST-2:L$$TUser@SQLSERVERDB",
        echo=False)
a_session = scoped_session(sessionmaker(autoflush=True, 
    bind=a_engine))
a_metadata = MetaData()
a_metadata.bind = a_engine
star = Table('stars', a_metadata, autoload=True)
opsim = Table('output_opsim3_61', a_metadata, autoload=True)
#galaxy = Table('galaxy', a_metadata, autoload=True)
#Galaxy = db.map(galaxy)
wd = Table('starsWD', a_metadata, autoload=True)
rrly = Table('starsRRLy', a_metadata, autoload=True)
bhb = Table('starsBHB', a_metadata, autoload=True)
db = SqlSoup(a_metadata)
Star = db.map(star)
OpSim3_61 = db.map(opsim)
Wd = db.map(wd)
BHB = db.map(bhb)
RRLy = db.map(rrly)

session = a_session
def initGalaxy(ra, dec, radiusdeg, columns, constraint=None):
    if constraint is not None:
        query = a_session.execute("EXECUTE [LSST].[dbo].[GalaxySearchSpecColsConstraint]\
            @RaSearch = %f, @DecSearch = %f, @apertureRadius = %f,\
            @ColumnNames = '%s', @WhereClause =\
            '%s'"%(ra,dec,radiusdeg*60.,columns,constraint))
    else:
        query = a_session.execute("EXECUTE [LSST].[dbo].[GalaxySearchSpecColsConstraint]\
            @RaSearch = %f, @DecSearch = %f, @apertureRadius = %f,\
            @ColumnNames = '%s'\
            '%s'"%(ra,dec,radiusdeg*60.,columns))
    coldesc = []
    for k in query.keys():
        coldesc.append({"name":k})
    return query, coldesc
'''
class Star(Entity):
  using_options(tablename="stars", autoload=True, metadata=a_metadata,
          session=a_session)
class Galaxy(Entity):
  using_options(tablename="galaxy", autoload=True, metadata=a_metadata,
          session=a_session)
class Wd(Entity):
  using_options(tablename="WDobs", autoload=True, metadata=a_metadata,
          session=a_session)

class OpSim3_61(Entity):
  using_options(tablename="output_opsim3_61", autoload=True,
          metadata=a_metadata, session=a_session)
'''
