#!/usr/bin/env python
from dbModel import *
from sqlalchemy import func
import datetime as dt
from pytz import timezone
import socket
import time

class LogEvents(object):
  def __init__(self, jobdescription="", jobid=None, ip = None):
    setup_all()
    create_all()
    self._tasknumber = None
    if jobid is None:
      jobid = b_session.query(func.max(CatalogEventLog.jobid)).one()[0]
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
    time.sleep(5)

  def persist(self, key, value, description):
    CatalogEventLog(jobid=self._jobid, pkey=unicode(key),
            pvalue=unicode(value),
            time=dt.datetime(1,1,1).now(timezone('US/Pacific')),
            taskNumber=self._tasknumber, ip=self._ip,
            description=unicode(description))
    b_session.commit()

  def registerTaskStart(self, tasknumber=None):
    key = "TASK_START" 
    if tasknumber is None:
      tasknumber = b_session.query(func.max(CatalogEventLog.taskNumber)).one()[0]
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
    setup_all()
    self._jobid = None
    self._states = {}
    if jobid is None:
      jobid = b_session.query(func.max(JobStateLog.jobid)).one()[0]
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
      statearr = JobStateLog.query.filter("jobid = %i and owner =\
              '%s'"%(self._jobid.getId(),self._jobid.getOwner())).all()
      for state in statearr:
        self._states[state.pkey] = state
    self.updateState("__Registration__", "Completed job registration")

  def getJobId(self):
    return self._jobid

  def getJobIdsByOwner(self, owner):
    jids = []
    idarr = b_session.query(func.distinct(JobStateLog.jobid)).filter("owner = '%s'"%(owner)).all()
    for id in idarr:
      jids.append(JobId(id[0], owner))
    return jids


  def updateState(self, key, state):
    if self._states.has_key(key):
      self._states[key].pvalue = unicode(state)
      self._states[key].time = dt.datetime(1,1,1).now(timezone('US/Pacific'))
      b_session.commit()
    else:
      self._states[key] = JobStateLog(jobid=self._jobid.getId(),
              owner=self._jobid.getOwner(), pkey=unicode(key),
              pvalue=unicode(state),
              time=dt.datetime(1,1,1).now(timezone('US/Pacific')))
      b_session.commit()

  def queryState(self, key):
    if self._states.has_key(key):
      b_session.refresh(self._states[key])
      return self._states[key].pvalue
    else:
      return None

  def showStates(self):
    states = {}
    keys = self._states.keys()
    for k in self._states.keys():
      b_session.refresh(self._states[k])
      #print k, self._states[k].pvalue
      states[k] = self._states[k].pvalue
    return states

  def refreshStates(self):
    for k in self._states.keys():
      b_session.refresh(self._states[k])

  def deleteStates(self):
    for key in self._states.keys():
      self._states[key].delete()
    b_session.commit()
