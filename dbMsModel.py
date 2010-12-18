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


session = scoped_session(sessionmaker())
a_engine =create_engine("mssql://LSST-2:L$$TUser@SQLSERVERDB",
        echo=False)
a_session = scoped_session(sessionmaker(autoflush=True, 
    bind=a_engine))
a_metadata = MetaData()
a_metadata.bind = a_engine
star = Table('stars', a_metadata, autoload=True)
opsim = Table('output_opsim3_61', a_metadata, autoload=True)
galaxy = Table('galaxy', a_metadata, autoload=True)
wd = Table('WDobs', a_metadata, autoload=True)
db = SqlSoup(a_metadata)
Star = db.map(star)
OpSim3_61 = db.map(opsim)
Galaxy = db.map(galaxy)
Wd = db.map(wd)

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
