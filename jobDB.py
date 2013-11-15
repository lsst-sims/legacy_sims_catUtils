#!/usr/bin/env python
from jobLogModel import engine, session, metadata,\
       CatalogEventLog, JobStateLog 
from sqlalchemy import func
import socket
import time
from datetime import tzinfo, timedelta, datetime


# A UTC class.

class UTC(tzinfo):
    _zero = timedelta(0)

    def utcoffset(self, dt):
                return self._zero

    def tzname(self, dt):
                return "UTC"

    def dst(self, dt):
                return self._zero


class LogEvents(object):
  def __init__(self, jobdescription="", jobid=None, ip = None):
    self._tasknumber = None
    self._conn = engine.connect()
    if jobid is None:
      jobid = session.query(func.max(CatalogEventLog.c.jobid)).one()[0]
      if jobid is None:
        self._jobid = 1
      else:
        self._jobid = jobid + 1
    else:
      self._jobid = int(jobid)
    self._jobdescription = jobdescription
    if ip is None:
      self._ip = socket.gethostbyname(socket.gethostname())
    else:
      self._ip = ip
    self.persist('__REGISTRATION__', 'Registered job %i'%self._jobid, '')

  def persist(self, key, value, description):
    ins = CatalogEventLog.insert().values(jobid=self._jobid, pkey=unicode(key),
            pvalue=unicode(value), time=datetime(1,1,1).now(UTC()),
            taskNumber=self._tasknumber, ip=self._ip,
            description=unicode(description))
    result = self._conn.execute(ins)

  def registerTaskStart(self, tasknumber=None):
    key = "TASK_START" 
    if tasknumber is None:
      tasknumber = session.query(func.max(CatalogEventLog.c.taskNumber)).one()[0]
      if tasknumber is None:
        tasknumber = 1
    else:
      pass
    self._tasknumber = tasknumber
    value = "Task started"
    self.persist(key, value, self._jobdescription)
    
  def registerEvent(self, eventid, eventdescription=""):
    key = "TASK_EVENT"
    value = eventid
    self.persist(key, value, eventdescription)

  def registerTaskStop(self, exitvalue=1):
    key = "TASK_STOP"
    value = "Task number %s stopped with term value %s"%(str(self._tasknumber), str(exitvalue))
    self.persist(key, value, self._jobdescription)

class JobId(object):
  def __init__(self, id, owner="anon"):
    self._jobid = int(id)
    self._owner = owner
  def setOwner(self, owner):
    self._owner = owner
  def getOwner(self):
    return self._owner
  def setId(self, id):
    self._jobid = id
  def getId(self):
    return self._jobid

class JobState(object):
  def __init__(self, jobid=None):
    self._jobid = None
    self._states = {}
    self._conn = engine.connect()
    if jobid is None:
      jobid = session.query(func.max(JobStateLog.c.jobid)).one()[0]
      if jobid is None:
        self._jobid = JobId(1)
      else:
        self._jobid = JobId(jobid + 1)
    else:
      if isinstance(jobid, JobId):
        self._jobid = jobid
      elif isinstance(jobid, int):
        self._jobid = JobId(jobid)
      else:
        raise Exception("The jobid is not an int or JobId class")
      statearr = session.query(JobStateLog).filter("jobid = %i and owner =\
              '%s'"%(self._jobid.getId(),self._jobid.getOwner())).all()
      for state in statearr:
        self._states[state.pkey] = state.pvalue
    self.updateState("__Registration__", "Completed job registration")

  def getJobId(self):
    return self._jobid

  def getJobIdsByOwner(self, owner):
    jids = []
    idarr = session.query(func.distinct(JobStateLog.c.jobid)).filter("owner = '%s'"%(owner)).all()
    for id in idarr:
      jids.append(JobId(id[0], owner))
    return jids


  def updateState(self, key, state):
    if self._states.has_key(key):
      self._states[key] = unicode(state)
      stmt = JobStateLog.update().where(
              JobStateLog.c.jobid == self._jobid.getId() and
              JobStateLog.c.owser == self._jobid.getOwner() and
              JobStateLog.c.pkey == unicode(key)).\
              values(pvalue=unicode(state), 
                     time=datetime.now(UTC()))
      result = self._conn.execute(stmt)  
    else:
      self._states[key] = unicode(state)
      ins = JobStateLog.insert().values(jobid=self._jobid.getId(), owner=self._jobid.getOwner(),
              pkey=unicode(key), pvalue=unicode(state), 
              time=datetime.now(UTC()))
      result = self._conn.execute(ins)

  def queryState(self, key):
    if self._states.has_key(key):
      return self._states[key]
    else:
      return None

  def showStates(self):
    states = {}
    for k in self._states:
      states[k] = self._states[k]
    return states

  def deleteStates(self):
    for key in self._states.keys():
      delstate = JobStateLog.delete().where(JobStateLog.c.jobid == self._jobid.getId() and
                                            JobStateLog.c.owner == self._jobid.getOwner() and 
                                            JobStateLog.c.pkey == unicode(key))
      result = self._conn.execute(delstate)
