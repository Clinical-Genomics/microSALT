""" Initial script to deliver and fetch data from a database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pymysql
import re
from microSALT import db
from microSALT.tables.samples import Samples
from microSALT.tables.seq_types import Seq_types
from microSALT.tables.profiles import Profiles
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
import sys
import yaml

import pdb # debug


#TODO: Rewrite all pushes/queries through session+commit
class DB_Manipulator:
 
  # Keeping mysql.yml seperate lets us share the main config one without security issues
  with open("{}/configs/mysql.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
    mysql = yaml.load(conf)

  def __init__(self, config, log):
    self.config = config
    self.logger = log
    self.engine = create_engine("mysql+pymysql://{}:{}@{}:{}/{}".format(self.mysql['user'],self.mysql['pw'],self.mysql['host'],self.mysql['port'],self.mysql['db']))
    Session = sessionmaker(bind=self.engine)
    self.session = Session()
    self.metadata = MetaData(self.engine)
    #TODO: Switch profiles to ORM format
    self.profiles = Profiles(self.metadata, self.config).tables

    self.create_tables()

  def create_tables(self):
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

  def get_100pc_alleles(self, cg_sid):
    """ Only grabs alleles from a samples that have 100% hit rate. Solves most problems"""
    return self.session.query(Seq_types.loci, Seq_types.allele).filter(Seq_types.CG_ID_sample==cg_sid, Seq_types.identity==100).all()

  def alleles2st(self, cg_sid):
    """Takes a CG sample ID and calculates most likely ST"""
    # TODO: Can be cleaned A LOT.
    organism = self.session.query(Samples.organism).filter(Samples.CG_ID_sample==cg_sid).scalar()
    #ONLY GRABS ALLELES WITH 100% ID
    alleles =self.get_100pc_alleles(cg_sid) 
    #Ugly check that duplicate entries dont exist
    alle_list = list()
    for item in alleles:
      if not item[0] in alle_list:
        alle_list.append(item[0])
      else:
        self.logger.error("Duplicate alleles found when converting alleles to ST. No logic implemented. Exited")
        sys.exit()
    for k,v in self.profiles.items():
      if k == organism:
        filtero = ""
        index = 0
        while index <  len(alleles):
          filtero += "v.c.{}=={}, ".format(alleles[index][0], alleles[index][1])
          index += 1
        ST = eval("self.session.query(v).filter({}).all()".format(filtero))
        if ST == []:
          self.logger.info("No ST for allele combo found. Setting ST to 0.")
          return 0
        if len(ST) > 1:
          self.logger.warning("Multiple ST found. Setting ST to -2. Manually investigate")
          return -2
        else: 
          return ST[0][0]
