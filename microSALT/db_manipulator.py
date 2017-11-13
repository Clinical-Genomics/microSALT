""" Initial script to deliver and fetch data from a database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pymysql
import re
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
import yaml

import pdb # debug

#TODO: Rewrite all pushes/queries through session+commit
class DB_Manipulator:

  # Keeping mysql.yml seperate lets us share the main config one without security issues
  with open("{}/mysql.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
    mysql = yaml.load(conf)

  def __init__(self, config, log):
    self.config = config
    self.logger = log
    self.engine = create_engine("mysql+pymysql://{}:{}@{}:{}/{}".format(self.mysql['user'],self.mysql['pw'],self.mysql['host'],self.mysql['port'],self.mysql['db']))
    Session = sessionmaker(bind=self.engine)
    self.session = Session()
    self.metadata = MetaData(self.engine)
    self.profiletables = dict()
    ## Tables
    self.projects = Table('projects', self.metadata,
        Column('CG_ID_project', String(15), primary_key=True),
        Column('date_analysis', DateTime),
        Column('organism', String(30)),
        Column('ST', SmallInteger, default=-1),
      )
    self.samples = Table('samples', self.metadata,
        Column('CG_ID_project', String(15), ForeignKey(self.projects.c.CG_ID_project), nullable=False),
        Column('CG_ID_sample', String(15), primary_key=True),
        Column('loci', String(10)),
        Column('assumed_ST', SmallInteger),
        Column('allele', SmallInteger),
        Column('haplotype', String(5)),
        Column('contig_name', String(20), primary_key=True),
        Column('contig_length', Integer),
        Column('contig_coverage', Float(6,2)),
        Column('identity', Float(3,2)),
        Column('evalue', String(10)),
        Column('bitscore', SmallInteger),
        Column('contig_start', Integer),
        Column('contig_end', Integer),
        Column('loci_start', Integer),
        Column('loci_end', Integer),
      )



    indata = os.listdir(self.config["folders"]["profiles"])
    for file in indata:
      #if not self.engine.dialect.has_table(self.engine, 'profile_{}'.format(file)):
        with open("{}/{}".format(self.config["folders"]["profiles"], file), "r") as fh:
          #Sets profile_* headers
          head = fh.readline()
          head = head.rstrip().split('\t')
          index = 0
          while index < len(head):
            if index == 0:
              header = "Table('profile_{}'.format(file), self.metadata,".format(file)
              header +="Column(head[{}], SmallInteger, primary_key=True),".format(index)
            else:
              header +="Column(head[{}], SmallInteger),".format(index)
            index = index+1
          header +=")"
          pt = eval(header)
          self.profiletables[file]=pt

    #Attempts to create tables every time a new object is formed. Maybe not ideal
    self.create_tables()

  def create_tables(self):
      if not self.engine.dialect.has_table(self.engine, 'projects'):
        self.projects.create()
        self.logger.info("Created projects table")
      if not self.engine.dialect.has_table(self.engine, 'samples'):
        self.samples.create()
        self.logger.info("Created samples table")
      for k,v in self.profiletables.items():
        if not self.engine.dialect.has_table(self.engine, "profile_{}".format(k)):
          self.profiletables[k].create()
          self.init_profiletable(k, v)       
 
  def add_record(self, data_dict, tablename):
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

  def update_record(self, data_dict, tablename, indict):
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
    return eval("dict.fromkeys(self.{}.c.keys())".format(table))

  def get_all_records(self, table):
    """ Returns all records for a given table"""
    return eval("self.{}.select().execute().fetchall()".format(table))


  def st2allele(self, organism, loci, ST):
    """ Takes organism name, loci name and assumed ST no.  and returns the correct allele number """
    for k,v in self.profiletables.items():
      if k == organism:
        try:
          return self.session.query(eval("v.c.{}".format(loci))).filter(v.c.ST==ST).scalar()
        except AttributeError:
          self.logger.info("Found loci {}, which has no profile entry".format(loci))
          return 0

  def alleles2st(self, cg_sid):
    """Takes a CG sample ID and calculates most likely ST"""
    # TODO: Can be cleaned A LOT.
    pids = self.session.query(self.samples.c.CG_ID_project).filter(self.samples.c.CG_ID_sample==cg_sid).all()
    if not pids.count(pids[0]) == len(pids):
      self.logger.error("One sample shared over multiple projects. Unable to set ST. Exited.")
      sys.exit()
    pid = pids[0][0]
    organism = self.session.query(self.projects.c.organism).filter(self.projects.c.CG_ID_project==pid).scalar()
    #ONLY GRABS ALLELES WITH 100% ID 
    alleles = self.session.query(self.samples.c.loci, self.samples.c.allele).filter(self.samples.c.CG_ID_sample==cg_sid, self.samples.c.identity==100).all()  
    #Ugly check that duplicate entries dont exist
    alle_list = list()
    for item in alleles:
      if not item[0] in alle_list:
        alle_list.append(item[0])
      else:
        self.logger.error("Duplicate alleles found when converting alleles to ST. No logic implemented. Exited")
        sys.exit()
    for k,v in self.profiletables.items():
      if k == organism:
        filtero = ""
        index = 0
        while index <  len(alleles):
          filtero += "v.c.{}=={}, ".format(alleles[index][0], alleles[index][1])
          index += 1
        ST = self.session.query(v).filter(exec(filtero)).all()
        if ST == []:
          self.logger.info("No ST for allele combo found. Setting ST to -1.")
          return 0
        else: 
          return ST


