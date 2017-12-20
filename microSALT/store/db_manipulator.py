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
    """Creates all tabes individually. A bit more control than usual"""
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

  def add_rec(self, data_dict, tablename):
    """Adds a record to the specified table through a dict with columns as keys."""
    table = eval(tablename)
    #Check for existing entry
    pk_list = table.__table__.primary_key.columns.keys()
    pk_values = list()
    for item in pk_list:
      pk_values.append(data_dict[item])
      existing = self.session.query(table).get(pk_values)

    #Add record
    if not existing:
      newobj = table()
      for k, v in data_dict.items():
        setattr(newobj, k, v)
      self.session.add(newobj)
      self.session.commit()
    else:
      self.logger.info("Failed to insert duplicate record into table {}".format(tablename))

  def upd_rec(self, data_dict, tablename, indict):
    """Updates a record to the specified table through a dict with columns as keys."""
    table = eval(tablename)
    filter = ""
    for k,v in data_dict.items():
      if v != None:
        filter +="{}=='{}' ,".format(k, v)
    if len(self.session.query(table).filter(filter).all()) > 1:
      self.logger.error("More than 1 record found when orm updating. Exited.")
      sys.exit()
    else:
      self.session.query(table).filter(filter).update(indict)
      self.session.commit()

  def init_profiletable(self, filename, table):
    """Creates profile tables by looping, since a lot of infiles exist"""
    data = table.insert()
    linedict = dict.fromkeys(table.c.keys())
    with open("{}/{}".format(self.config["folders"]["profiles"], filename), "r") as fh:
      #Skips header
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
  
  def get_columns(self, tablename):
    """ Returns all records for a given ORM table"""
    table = eval(tablename)
    return dict.fromkeys(table.__table__.columns.keys())

  def setPredictor(self, cg_sid, pks=dict()):
    """ Helper function that flags a set of seq_types as part of the final prediction.
      Flags all in case nothing can be established.
      Optionally takes in a loci -> contig name dict for allele distinction. """
    sample = self.session.query(Seq_types).filter(Seq_types.CG_ID_sample==cg_sid)
    if pks == dict():
      sample.update({Seq_types.st_predictor: 1})
    else:
      for loci, cn in pks.items():
        sample.filter(Seq_types.loci==loci, Seq_types.contig_name==cn).update({Seq_types.st_predictor : 1})
    self.session.commit()

  def alleles2st(self, cg_sid):
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    if organism is None:
      self.logger.warning("No organism set for {}. Most likely control sample. Setting ST to -1".format(cg_sid))
      return -1
    hits = self.session.query(Seq_types.loci, Seq_types.allele, Seq_types.identity).filter(Seq_types.CG_ID_sample==cg_sid, Seq_types.identity>=99.9, Seq_types.evalue==0.0).all()
    #Unused, might be valuable later
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
        self.setPredictor(cg_sid)
        return -3

    # Tests all allele combinations found to see if any of them result in ST
    filter = ""
    for key, val in uniqueDict.items():
      subfilter = ""
      for num in val:
        subfilter += " self.profiles[organism].c.{}=={} ,".format(key, num)
      if len(num) > 1:
        subfilter = "_or({})".format(subfilter)
      filter += subfilter
    filter = "_and({})".format(filter)

    output = self.session.query(self.profiles[organism]).filter(text(filter)).all()
    if len(output) > 1:
      STlist= list()
      for st in output:
        STlist.append(st.ST)
      best = self.bestST(cg_sid, STlist)
      self.logger.warning("Multiple possible ST found for sample {}, list: {}. Established ST{} as best hit.".format(cg_sid, STlist, best))
      return best
    elif len(output) == 1:
      #This bestST call is excessive and should be streamlined later
      return self.bestST(cg_sid, [output[0].ST])
    else:
      self.logger.warning("Sample {} on {} has a single allele set but no matching ST. Either incorrectly called allele, or novel ST has been found. Setting ST to -2".format(cg_sid, organism))
      self.setPredictor(cg_sid)
      return -2

  def bestST(self, cg_sid, st_list):
    """Takes in a list of ST and a sample. Establishes which ST is most likely by criteria  id -> eval -> contig coverage"""

    profiles = list()
    scores = dict()
    bestalleles = dict()
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    for st in st_list:
      scores[st] = dict()
      bestalleles[st] = dict()
      scores[st]['id'] = 0
      scores[st]['eval'] = 0
      scores[st]['cc'] = 0
      profiles.append(self.session.query(self.profiles[organism]).filter(text('ST={}'.format(st))).first())

    # Get value metrics for each allele set that resolves an ST
    # TODO: Make it readable 
    for prof in profiles:
      filterstring = "and_(Seq_types.CG_ID_sample=='{}', or_(".format(cg_sid)
      for index, allele in enumerate(prof):
        if 'ST' not in prof.keys()[index] and 'clonal_complex' not in prof.keys()[index]:
          filterstring += "and_(Seq_types.loci=='{}', Seq_types.allele=='{}'), ".format(prof.keys()[index], allele)
      filterstring = filterstring[:-2] + "))"
      alleles = self.session.query(Seq_types).filter(eval(filterstring)).all()
   

      # Keep only best allele of each loci name
      index = 0
      seenLoci = dict()
      while index < len(alleles)-1:
        if alleles[index].loci in seenLoci:
          prev = seenLoci[alleles[index].loci]

          if (alleles[index].identity > alleles[prev].identity) or\
          (alleles[index].identity == alleles[prev].identity and float(alleles[index].evalue) < float(alleles[prev].evalue))\
          or (alleles[index].identity == alleles[prev].identity and float(alleles[index].evalue) == float(alleles[prev].evalue)\
          and alleles[index].contig_coverage > alleles[prev].contig_coverage):
            alleles.pop(prev)
            seenLoci[alleles[index].loci] = index
          else:
            alleles.pop(index)
        else:    
          seenLoci[alleles[index].loci] = index
          index += 1

      #Compute metrics
      for allele in alleles:
        scores[prof.ST]['id'] += allele.identity
        scores[prof.ST]['eval'] += float(allele.evalue)
        scores[prof.ST]['cc'] += allele.contig_coverage
        bestalleles[prof.ST][allele.loci] = ""
        bestalleles[prof.ST][allele.loci] += allele.contig_name

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
    self.setPredictor(cg_sid, bestalleles[topST])
    return topST
