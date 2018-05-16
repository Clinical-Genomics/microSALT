""" Delivers and fetches data from the database
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import sys
import warnings

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker

from microSALT.store.orm_models import app, Projects, Resistances, Samples, Seq_types, Versions
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
    #Turns off pymysql deprecation warnings until they can update their code
    with warnings.catch_warnings():
      warnings.simplefilter("ignore")
      self.create_tables()

  def create_tables(self):
    """Creates all tables individually. A bit more control than usual"""
    if not self.engine.dialect.has_table(self.engine, 'projects'):
      Projects.__table__.create(self.engine)
      self.logger.info("Created projects table")
    if not self.engine.dialect.has_table(self.engine, 'samples'):
      Samples.__table__.create(self.engine)
      self.logger.info("Created samples table")
    if not self.engine.dialect.has_table(self.engine, 'versions'):
      Versions.__table__.create(self.engine)
      self.logger.info("Created versions table")
    if not self.engine.dialect.has_table(self.engine, 'seq_types'):
      Seq_types.__table__.create(self.engine)
      self.logger.info("Created sequencing types table")
    if not self.engine.dialect.has_table(self.engine, 'resistances'):
      Resistances.__table__.create(self.engine)
      self.logger.info("Created resistance table")
    for k,v in self.profiles.items():
      if not self.engine.dialect.has_table(self.engine, "profile_{}".format(k)):
        self.profiles[k].create()
        self.init_profiletable(k, v)       
        self.add_rec({'name': 'profile_{}'.format(k), 'version': '0'}, 'Versions', force=True)
        self.logger.info("Profile table profile_{} initialized".format(k))

  def add_rec(self, data_dict, tablename, force=False):
    """Adds a record to the specified table through a dict with columns as keys."""
    try:
      table = eval(tablename)
      #Check for existing entry
      pk_list = table.__table__.primary_key.columns.keys()
    except Exception as e:
      self.logger.error("Attempted to access table {} which has not been created".format(tablename))
    pk_values = list()
    for item in pk_list:
      pk_values.append(data_dict[item])
    existing = self.session.query(table).get(pk_values)
    #Add record
    if not existing or force:
      newobj = table()
      for k, v in data_dict.items():
        setattr(newobj, k, v)
      self.session.add(newobj)
      self.session.commit()
    else:
      self.logger.warn("Hindered insertion of existing record ({}) into table {}".format(data_dict, tablename))

  def upd_rec(self, req_dict, tablename, upd_dict):
    """Updates a record to the specified table through a dict with columns as keys."""
    table = eval(tablename)
    args = list()
    for k,v in req_dict.items():
      if v != None:
        args.append("table.{}=='{}'".format(k, v))
    filter = ','.join(args)
    if len(self.session.query(table).filter(eval(filter)).all()) > 1:
      self.logger.error("More than 1 record found when orm updating. Exited.")
      sys.exit()
    else:
      self.session.query(table).filter(eval(filter)).update(upd_dict)
      self.session.commit()

  def reload_profiletable(self, organism):
    """Drop the named non-orm table, then load it with fresh data"""
    table = self.profiles[organism]
    self.profiles[organism].drop()
    self.profiles[organism].create()
    self.init_profiletable(organism, table)
    self.logger.info("Profile table for {} updated to latest version".format(organism)) 
 
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
 
  def get_columns(self, tablename):
    """ Returns all records for a given ORM table"""
    table = eval(tablename)
    return dict.fromkeys(table.__table__.columns.keys())

  def get_profiles(self):
    """Returns non-orm profiles tables"""
    return self.profiles

  def exists(self, table, item):
    """ Takes a k-v pair and checks for the entrys existence in the given table """
    filterstring = ""
    for k, v in item.items():
      filterstring += "{}.{}=='{}',".format(table, k, v)
    filterstring = filterstring[:-1]
    table = eval(table)
    entry = self.session.query(table).filter(eval(filterstring)).scalar()
    if entry is None:
      return False
    else:
      return True

  def get_version(self, name):
    """ Gets the version from a given name. Should be generalized to return any value for any input"""
    version = self.session.query(Versions).filter(Versions.name==name).scalar()
    if version is None:
      return "0"
    else:
      return version.version

  def setPredictor(self, cg_sid, pks=dict()):
    """ Helper function. Flags a set of seq_types as part of the final prediction.
    Uses optional pks[loci]=contig_name dictionary to distinguish in scenarios where an allele number has multiple hits"""
    sample = self.session.query(Seq_types).filter(Seq_types.CG_ID_sample==cg_sid)

    if pks == dict():
      sample.update({Seq_types.st_predictor: 1})
    else:
      #Resets all previous predictors
      sample.update({Seq_types.st_predictor: None})
      #Set subset
      for loci, cn in pks.items():
        if isinstance(cn, list):
          for entry in cn:
            sample.filter(Seq_types.loci==loci, Seq_types.allele==entry).update({Seq_types.st_predictor : 1})
        else:
          sample.filter(Seq_types.loci==loci, Seq_types.contig_name==cn).update({Seq_types.st_predictor : 1})
    self.session.commit()

  def alleles2st(self, cg_sid):
    """ Takes a CG_ID_sample and predicts the correct ST """
    threshold = True
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    if organism is None:
      self.logger.warning("No organism set for {}. Most likely control sample. Setting ST to -1".format(cg_sid))
      return -1
    [alleles, allelediff] = self.get_unique_alleles(cg_sid, organism, threshold)
    if allelediff < 0:
      threshold = False
      [alleles, allelediff] = self.get_unique_alleles(cg_sid, organism, threshold)
      if allelediff < 0:
        self.logger.warning("Insufficient allele hits to establish ST for sample {}, even without thresholds. Setting ST to -3"\
                            .format(cg_sid, organism))
        self.setPredictor(cg_sid)
        return -3
    self.upd_rec({'CG_ID_sample':cg_sid}, 'Samples', {'aux_alleles':allelediff})

    # Tests all allele combinations found to see if any of them result in ST
    filter = list()
    for key, val in alleles.items():
      subfilter = list()
      for num in val:
        subfilter.append(" self.profiles[organism].c.{}=={} ".format(key, num))
      subfilter = ','.join(subfilter)
      if len(val) > 1:
        subfilter = "or_({})".format(subfilter)
      filter.append(subfilter)
    filter = ','.join(filter)
    filter = "and_({})".format(filter)
    output = self.session.query(self.profiles[organism]).filter(eval(filter)).all()
 
    if len(output) > 1:
      STlist= list()
      for st in output:
        STlist.append(st.ST)
      best = self.bestST(cg_sid, STlist)
      self.upd_rec({'CG_ID_sample':cg_sid}, 'Samples', {'aux_ST':1})
      if threshold:   
        self.logger.warning("Multiple ST within threshold found for sample {}, list: {}. Established ST{} as best hit.".format(cg_sid, STlist, best))
      return best
    elif len(output) == 1:
      #Doing bestST only to establish best loci number combination
      return self.bestST(cg_sid, [output[0].ST])
    elif threshold:
      self.logger.warning("Sample {} on {} has an allele set but no ST. Novel ST found, setting ST to -4".format(cg_sid, organism))
      bestSet = self.bestAlleles(cg_sid)
      self.setPredictor(cg_sid, bestSet)
      return -4
    else:
      self.logger.warning("Sample {} on {} has an allele set but hits are low-quality and\
 do not resolve to an ST. Setting ST to -2".format(cg_sid, organism))
      bestSet = self.bestAlleles(cg_sid)
      self.setPredictor(cg_sid, bestSet)
      return -2

  def bestST(self, cg_sid, st_list):
    """Takes in a list of ST and a sample. Establishes which ST is most likely by criteria id*span -> eval -> contig coverage"""
    profiles = list()
    scores = dict()
    bestalleles = dict()
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    for st in st_list:
      scores[st] = dict()
      bestalleles[st] = dict()
      scores[st]['spanid'] = 0
      scores[st]['eval'] = 0
      scores[st]['cc'] = 0
      scores[st]['span'] = 0
      profiles.append(self.session.query(self.profiles[organism]).filter(text('ST={}'.format(st))).first())

    # Get values for each allele set that resolves an ST
    for prof in profiles:
      alleleconditions = list()
      alleledict = dict()
      allconditions = ["Seq_types.CG_ID_sample=='{}'".format(cg_sid)]

      for index, allele in enumerate(prof):
        if 'ST' not in prof.keys()[index] and 'clonal_complex' not in prof.keys()[index]:
          condition = "Seq_types.loci=='{}' , Seq_types.allele=='{}'".format(prof.keys()[index], allele)
          alleledict[prof.keys()[index]] = ""
          alleleconditions.append("and_({})".format(condition))

      alleleconditions = "or_({})".format(','.join(alleleconditions))
      allconditions.append(alleleconditions)
      allconditions = "and_({})".format(','.join(allconditions))
      all_alleles = self.session.query(Seq_types).filter(eval(allconditions)).all()

      # Keep only best hit each loci
      for allele in all_alleles:
       if alleledict[allele.loci] == "":
         alleledict[allele.loci] = allele
       else:
        if ((allele.span*allele.identity > alleledict[allele.loci].span*alleledict[allele.loci].identity) or\
        (allele.span*allele.identity == alleledict[allele.loci].span*alleledict[allele.loci].identity and\
        float(allele.evalue) < float(alleledict[allele.loci].evalue)) or\
        (allele.span*allele.identity == alleledict[allele.loci].span*alleledict[allele.loci].identity and\
        float(allele.evalue) == float(alleledict[allele.loci].evalue) and\
        allele.contig_coverage > alleledict[allele.loci].contig_coverage)):
          alleledict[allele.loci] = allele

      #Create score dict for the ST
      for key, allele in alleledict.items():
        scores[prof.ST]['spanid'] += allele.span*allele.identity
        scores[prof.ST]['eval'] += float(allele.evalue)
        scores[prof.ST]['cc'] += allele.contig_coverage
        bestalleles[prof.ST][allele.loci] = ""
        bestalleles[prof.ST][allele.loci] += allele.contig_name

    #Establish best ST
    topST = ""
    topID = 0
    topEval = 100
    topCC = 0
    for key, val in scores.items():
      if scores[key]['spanid'] > topID:
        topID = scores[key]['spanid']
        topST = key
      elif scores[key]['spanid'] == topID and scores[key]['eval'] < topEval:
        topEval = scores[key]['eval']
        topST = key
      elif scores[key]['spanid'] == topID and scores[key]['eval'] == topEval:
        topCC = scores[key]['cc']
        topST = key
    self.setPredictor(cg_sid, bestalleles[topST])
    return topST

  def bestAlleles(self, cg_sid):
    """ Establishes which allele set (for novel ST) is most likely by criteria span* id -> eval -> contig coverage"""
    hits = self.session.query(Seq_types.contig_name, Seq_types.loci, Seq_types.span, Seq_types.identity, Seq_types.evalue, Seq_types.contig_coverage)\
           .filter(Seq_types.CG_ID_sample==cg_sid).all()
    bestHits = dict()
    alleledict = dict()
    for allele in hits:
      if allele.loci not in bestHits.keys():
        bestHits[allele.loci] = allele.contig_name
        alleledict[allele.loci] = [allele.identity, allele.evalue, allele.contig_coverage, allele.span]
      else:
        if ((allele.identity*allele.span > alleledict[allele.loci][0]*alleledict[allele.loci][3]) or\
        (allele.identity*allele.span == alleledict[allele.loci][0]*alleledict[allele.loci][3] and\
        float(allele.evalue) < float(alleledict[allele.loci][1])) or\
        (allele.identity*allele.span == alleledict[allele.loci][0]*alleledict[allele.loci][3] and\
        float(allele.evalue) == float(alleledict[allele.loci][1]) and\
        allele.contig_coverage > alleledict[allele.loci][2])):
          bestHits[allele.loci] = allele.contig_name
          alleledict[allele.loci] = [allele.identity, allele.evalue, allele.contig_coverage, allele.span]
    return bestHits


  def get_unique_alleles(self, cg_sid, organism, threshold=True):
    """ Returns a dict containing all unique alleles at every loci, and allele difference from expected"""
    if threshold:
      hits = self.session.query(Seq_types.loci, Seq_types.allele)\
           .filter(Seq_types.CG_ID_sample==cg_sid).all()
    else:
      hits = self.session.query(Seq_types.loci, Seq_types.allele).filter(Seq_types.CG_ID_sample==cg_sid).all()

    #Establish number of unique hits
    uniqueDict = dict()
    for hit in hits:
      if hit.loci not in uniqueDict.keys():
        uniqueDict[hit.loci] = list()
        uniqueDict[hit.loci].append(hit.allele)
      elif hit.allele not in uniqueDict[hit.loci]:
        uniqueDict[hit.loci].append(hit.allele)

    if 'clonal_complex' in self.profiles[organism].columns.keys():
      non_allele_columns = 2
    else:
      non_allele_columns = 1
    allele_overabundance = len(uniqueDict.keys()) - (len(self.profiles[organism].columns.values()) - non_allele_columns)
    return [uniqueDict, allele_overabundance]
