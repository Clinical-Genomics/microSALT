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

  def alleles2st(self, cg_sid):
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    hits = self.session.query(Seq_types.loci, Seq_types.allele, Seq_types.identity).filter(Seq_types.CG_ID_sample==cg_sid, Seq_types.identity>=99.9, Seq_types.evalue==0.0).all()
    thresholdless = True
    if 'clonal_complex' in self.profiles[organism].columns.keys():
      non_allele_columns = 2
    else:
      non_allele_columns = 1 

    #Establish there's enough unique hits at all, otherwise we're screwed from the get go
    uniqueDict = dict()
    for hit in hits:
      if hit.loci not in uniqueDict.keys():
        uniqueDict[hit.loci] = list()
        uniqueDict[hit.loci].append(hit.allele)
      elif hit.allele not in uniqueDict[hit.loci]:
        uniqueDict[hit.loci].append(hit.allele)
    if len(self.profiles[organism].columns.values()) - non_allele_columns > len(uniqueDict.keys()):
      # Not enough hits with thresholds, using threshold-less. 
      hits = self.session.query(Seq_types.loci, Seq_types.allele, Seq_types.identity).filter(Seq_types.CG_ID_sample==cg_sid).all()
      thresholds = False
      #If Still not enough, return -3
      for hit in hits:
        if hit.loci not in uniqueDict.keys():
          uniqueDict[hit.loci] = list()
          uniqueDict[hit.loci].append(hit.allele)
        elif hit.allele not in uniqueDict[hit.loci]:
          uniqueDict[hit.loci].append(hit.allele)
      if len(self.profiles[organism].columns.values()) - non_allele_columns > len(uniqueDict.keys()):
        self.logger.warning("Insufficient allele hits to establish ST for sample {}, even without thresholds. Setting ST to -3".format(cg_sid, organism))
        return -3

    # Tests all allele combinations found to see if any of them result in ST
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
      best = self.bestST(STlist, cg_sid)
      self.logger.warning("Multiple possible ST found for sample {}, list: {}. Established ST{} as best hit.".format(cg_sid, STlist, best))
      return best
    elif len(output) == 1:
      return output[0].ST
    else:
      self.logger.warning("Sample {} on {} has a single allele set but no matching ST. Either incorrectly called allele, or novel ST has been found. Setting ST to -2".format(cg_sid, organism))
      return -2

  def bestST(self, st_list, cg_sid):
    """Takes in a list of ST and a sample. Establishes which ST is most likely by criteria  id -> eval -> contig coverage"""

    profiles = list()
    scores = dict()
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    for st in st_list:
      scores[st] = dict()
      scores[st]['id'] = 0
      scores[st]['eval'] = 0
      scores[st]['cc'] = 0
      profiles.append(self.session.query(self.profiles[organism]).filter(text('ST={}'.format(st))).first())

    # Get value metrics for each allele set that resolves an ST 
    for prof in profiles:
      filterstring = "and_(Seq_types.CG_ID_sample=='{}', or_(".format(cg_sid)
      for index, allele in enumerate(prof):
        if 'ST' not in prof.keys()[index] and 'clonal_complex' not in prof.keys()[index]:
          filterstring += "and_(Seq_types.loci=='{}', Seq_types.allele=='{}'), ".format(prof.keys()[index], allele)
      filterstring = filterstring[:-2] + "))"
      alleles = self.session.query(Seq_types).filter(eval(filterstring)).all()

      # Keep only best allele of each loci name
      index = 0
      addedloci = dict()
      while index < len(alleles)-1:
        if alleles[index].loci in addedloci:
          nex = addedloci[alleles[index].loci]

          if (alleles[index].identity > alleles[nex].identity) or (alleles[index].identity == alleles[nex].identity and float(alleles[index].evalue) < float(alleles[nex].evalue))\
          or (alleles[index].identity == alleles[nex].identity and float(alleles[index].evalue) == float(alleles[nex].evalue) and alleles[index].contig_coverage > alleles[nex].contig_coverage):
            alleles.pop(nex)
            addedloci[alleles[index].loci] = index
          else:
            alleles.pop(index)
        else:    
          addedloci[alleles[index].loci] = index
          index += 1
      #Compute metrics
      for allele in alleles:
        scores[prof.ST]['id'] += allele.identity
        scores[prof.ST]['eval'] += float(allele.evalue)
        scores[prof.ST]['cc'] += allele.contig_coverage

    topST = ""
    topID = 0
    topEval = 100
    topCC = 0
    for key, val in scores.items():
      if scores[key]['id'] > topID:
        topID = scores[key]['id']
        topST = key
      elif scores[key]['id'] == topID and scores[key]['eval'] < topEval:
        topEval = scores[key]['eval']
        topST = key
      elif scores[key]['id'] == topID and scores[key]['eval'] == topEval:
        topCC = scores[key]['cc']
        topST = key
    return topST
