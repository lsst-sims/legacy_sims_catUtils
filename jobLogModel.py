import warnings
from sqlalchemy.pool import NullPool
from sqlalchemy.orm import scoped_session, sessionmaker
from sqlalchemy.sql import expression
from sqlalchemy import create_engine, Column, Table, Integer,\
                       UnicodeText, DateTime,\
                       Text
from sqlalchemy import ThreadLocalMetaData
import sqlalchemy.databases as sd 
from sqlalchemy import func
from sqlalchemy import schema
from sqlalchemy import exc as sa_exc

warnings.simplefilter("ignore", category=sa_exc.SAWarning)
engine = create_engine("sqlite:///tmp.sqlite",
        echo=False, convert_unicode=False, poolclass=NullPool)
session = application_session = scoped_session(sessionmaker(autoflush=True,
     bind=engine))
metadata = ThreadLocalMetaData()
metadata.bind = engine

CatalogEventLog = Table('eventlog', metadata,
            Column('jobid', Integer, index=True),
            Column('owner', UnicodeText),
            Column('pkey', UnicodeText),
            Column('pvalue', UnicodeText),
            Column('time', DateTime(timezone=True)),
            Column('taskNumber', Integer),
            Column('ip', Text),
            Column('description', UnicodeText))
try:
    CatalogEventLog.create(engine)
except sa_exc.OperationalError:
    #Catalog already exists
    pass

JobStateLog = Table('statelog', metadata,
        Column('jobid', Integer, index=True),
        Column('owner', Text),
        Column('pkey', Text),
        Column('pvalue', Text),
        Column('time', DateTime(timezone=True)))
try:
    JobStateLog.create(engine)
except sa_exc.OperationalError:
    #Database already exists
    pass
