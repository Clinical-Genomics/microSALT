""" Initial script to deliver and fetch data from a database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import re
import sys
import yaml

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker

from microSALT import app
from microSALT.store.orm_models import Projects, Samples, Seq_types
from microSALT.store.models import Profiles
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker

#TODO: Rewrite all pushes/queries through session+commit
#TODO: Contains a lot of legacy from deving. Remove what cant be repurposed.
#TODO: Direct most server calls through here.
class DB_Manipulator:
 
  def __init__(self, config, log):
    self.config = config
    self.logger = log
    self.engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'])
    Session = sessionmaker(bind=self.engine)
    self.session = Session()
    self.metadata = MetaData(self.engine)
    #TODO: Switch profiles to ORM format
    self.profiles = Profiles(self.metadata, self.config).tables

    self.create_tables()

  def create_tables(self):
      if not self.engine.dialect.has_table(self.engine, 'projects'):
        Projects.__table__.create(self.engine)
        self.logger.info("Created samples table")
      if not self.engine.dialect.has_table(self.engine, 'samples'):
        Samples.__table__.create(self.engine)
        self.logger.info("Created samples table")
      if not self.engine.dialect.has_table(self.engine, 'seq_types'):
        Seq_types.__table__.create(self.engine)
        self.logger.info("Created sequencing types table")
      for k,v in self.profiles.items():
        if not self.engine.dialect.has_table(self.engine, "profile_{}".format(k)):
          self.profiles[k].create()
          self.init_profiletable(k, v)       

  def _is_orm(self, tablename):
    try:
      eval("self.{}.__table__".format(tablename))
    except AttributeError:
      return False
 
  def add_rec(self, data_dict, tablename):
    """ Adds a record to the specified table through a dict with columns as keys.
        Use in tandem with get_columns for best results """
   
    #Find PK 
    pk_list = [key.name for key in inspect(eval("self.{}".format(tablename))).primary_key]
    #Check key collision, REALLY have to streamline this.
    statem  = "self.session.query("
    for item in pk_list:
      statem += "self.{}.c.{}, ".format(tablename, item)
    statem +=").filter("
    for item in pk_list:
      statem += "self.{}.c.{}=='{}', ".format(tablename, item, data_dict[item])
    statem +=").all()"

    if not eval(statem):
      eval("self.{}.insert().execute(data_dict)".format(tablename))
      self.session.commit() 
    else:
      self.logger.info("Failed to insert duplicate record into table {}".format(tablename))

  def add_rec_orm(self, data_dict, tablename):
    #Check for duplicate
    pk_list = eval("{}.__table__.primary_key.columns.keys()".format(tablename))
    #Check key collision, REALLY have to streamline this.
    statem  = "self.session.query("
    for item in pk_list:
      statem += "{}.{}, ".format(tablename, item)
    statem +=").filter("
    for item in pk_list:
      statem += "{}.{}=='{}', ".format(tablename, item, data_dict[item])
    statem +=").all()"

    if not eval(statem):   
      newobj = eval("{}()".format(tablename))
      for k, v in data_dict.items():
        exec("newobj.{} = v".format(k, v))
      self.session.add(newobj)
      self.session.commit()
    else:
      self.logger.info("Failed to insert duplicate record into table {}".format(tablename))

  def upd_rec(self, data_dict, tablename, indict):
    #Oh dear.. Just happy it works tbh
    instring = ""
    for k,v in data_dict.items():
      #Avoids invalid match with empty fields
      if not v==None:
        instring+="self.{}.c.{}=='{}', ".format(tablename, k,v) 
    if len(eval("self.session.query(self.{}).filter({}).all()".format(tablename, instring))) > 1:
      self.logger.error("More than 1 record found when updating. Exited.")
      sys.exit()
    input = "self.{}.update(self.{}).where(".format(tablename,tablename)
    for k,v in data_dict.items():
      input +="self.{}.c.{}=='{}' and ".format(tablename, k,v)
    input = input[:-4] + ").values("
    for k,v in indict.items(): 
      input += "{}={}".format(k,v)
    input += ")"
    eval(input).execute()

  def upd_rec_orm(self, data_dict, tablename, indict):
    argstring = ""
    for k,v in data_dict.items():
      if v != None:
        argstring +="{}='{}', ".format(k, v)
    if len(eval("self.session.query({}).filter_by({}).all()".format(tablename, argstring))) > 1:
      self.logger.error("More than 1 record found when orm updating. Exited.")
      sys.exit()
    else:
      eval("self.session.query({}).filter_by({}).update({})".format(tablename, argstring, indict))  
      self.session.commit()

  def init_profiletable(self, filename, table):
    data = table.insert()
    linedict = dict.fromkeys(table.c.keys())
    with open("{}/{}".format(self.config["folders"]["profiles"], filename), "r") as fh:
      #Skip header
      head = fh.readline()
      head = head.rstrip().split('\t')
      for line in fh:
        line = line.rstrip().split('\t')
        index = 0
        while index < len(line):
          linedict[head[index]] = line[index]
          index = index+1
        data.execute(linedict)
    self.logger.info("Created profile table {}".format(table))
  
  def get_columns(self, table):
    """ Returns a dictionary where each column is a key. Makes re-entry easier """
    return eval("dict.fromkeys(self.{}.keys())".format(table))

  def get_columns_orm(self, table):
    """ Returns all records for a given ORM table"""
    return eval("dict.fromkeys({}.__table__.columns.keys())".format(table))

  def get_all_records(self, table):
    """ Returns all records for a given table"""
    return eval("self.{}.select().execute().fetchall()".format(table))


  def st2allele(self, organism, loci, ST):
    """ Takes organism name, loci name and assumed ST no.  and returns the correct allele number """
    #LEGACY CODE. Keeping for the moment in case need re-emerges.
    for k,v in self.profiles.items():
      if k == organism:
        try:
          return self.session.query(eval("v.c.{}".format(loci))).filter(v.c.ST==ST).scalar()
        except AttributeError:
          self.logger.info("Found loci {}, which has no profile entry".format(loci))
          return 0

  def simpleAllele2st(self, cg_sid):
    """ Finds a definitive ST if such exists """
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    # Find clear top hits
    hits = self.session.query(Seq_types.loci, Seq_types.allele, Seq_types.identity).filter(Seq_types.CG_ID_sample==cg_sid, Seq_types.identity>=99.9, Seq_types.evalue==0.0).all()

    #Check that there's only correct number of loci hits
    if len(self.profiles[organism].columns.values()) == len(hits):
      self.logger.info("Equal number of hits as organism profile for {} on {}".format(cg_sid, organism))

    #Find ST
    filterstring = ""
    for entry in hits:
      filterstring += " {}={} and".format(entry.loci, entry.allele)
    filterstring = filterstring[:-3]
    output = self.session.query(self.profiles[organism]).filter(text(filterstring)).all()

    if len(output) == 1:
      return output[0].ST
    else:
      return -1 #No Simple solution found

  def complexAllele2st(self, cg_sid):
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    if 'clonal_complex' in self.profiles[organism].columns.keys():
      non_allele_columns = 2
    else:
      non_allele_columns = 1 
    #All hits
    hits = self.session.query(Seq_types.loci, Seq_types.allele, Seq_types.identity).filter(Seq_types.CG_ID_sample==cg_sid).all()

    #Establish there's enough unique hits at all, otherwise we're screwed from the get go
    uniqueDict = dict()
    for hit in hits:
      if hit.loci not in uniqueDict.keys():
        uniqueDict[hit.loci] = list()
        uniqueDict[hit.loci].append(hit.allele)
      elif hit.allele not in uniqueDict[hit.loci]:
        uniqueDict[hit.loci].append(hit.allele)

    #Check that correct amount of alleles has been found AND that no allele has multiple numbers
    #This corrects the issue where the loci does not have 100% id
    if len(self.profiles[organism].columns.values()) - non_allele_columns == len(uniqueDict.keys()) and sum([len(x) for x in uniqueDict.values()]) == len(self.profiles[organism].columns.values()) - non_allele_columns:
      # Find ST (extended)
      filterstring = ""
      for entry in hits:
        filterstring += " {}={} and".format(entry.loci, entry.allele)
      filterstring = filterstring[:-3]
      output = self.session.query(self.profiles[organism]).filter(text(filterstring)).all()
  
      if len(output) > 0:
        self.logger.warning("Sampe {} has ST {} from extended search (breaking thresholds). Be aware".format(cg_sid, output[0].ST))
        return output[0].ST
      else:
        self.logger.warning("Sufficent loci found but no ST. True hit is not top hit, or novel ST has been found. Sample {} on {}. Setting ST to -2".format(cg_sid, organism))
        return -2
    elif len(self.profiles[organism].columns.values()) - non_allele_columns > len(uniqueDict.keys()):
      self.logger.warning("Less hits than organism profile for {} on {}. Relax blast constraint? Unable to fanthom a ST. Setting to -3".format(cg_sid, organism))
      return -3

    # This corrects the issue where multiple alleles for a given loci were found.
    # Find all ST combinations, put them in a list
    filterstring = ""
    for key, val in uniqueDict.items():
      for num in val:
        if val.index(num) == 0 and len(val) > 1:
          filterstring += " ("

        filterstring += " {}={} ".format(key, num) 

        if len(val) == 1:
          filterstring += "and"
        elif val.index(num)+1 == len(val) and len(val) > 1:
          filterstring += ") and"
        else:
          filterstring += "or"
    filterstring = filterstring[:-3]
    output = self.session.query(self.profiles[organism]).filter(text(filterstring)).all()
    if len(output) > 1:
      STlist= list()
      for st in output:
        STlist.append(st.ST) 
      self.logger.warning("Multiple possible ST found for sample {}, list: {}".format(cg_sid, STlist))
    else:
      self.logger.warning("Established ST{} for {} from expandeded search".format(output[0].ST, cg_sid))
      return output[0].ST

  def alleles2st(self, cg_sid):
    """Takes a CG sample ID and calculates most likely ST"""
    # All this needs to be conveyed to views
    ST = self.simpleAllele2st(cg_sid)
    if ST > 0:
      return ST
    #TODO: Distinquish between different negative codes
    else:
      self.logger.warning("No Simple ST found for {}. Using complex finder.".format(cg_sid))
      return self.complexAllele2st(cg_sid)
