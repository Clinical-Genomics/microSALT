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

  def __init__(self, config):
    self.config = config
    self.engine = create_engine("mysql+pymysql://{}:{}@{}:{}/{}".format(self.mysql['user'],self.mysql['pw'],self.mysql['host'],self.mysql['port'],self.mysql['db']))
    Session = sessionmaker(bind=self.engine)
    self.session = Session()
    self.metadata = MetaData(self.engine)
    self.profiletables = dict()
    ## Tables
    self.samples = Table('samples', self.metadata,
        Column('CG_ID_Sample', String(15), primary_key=True),
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
    #Assuming this ID mess holds, still unsure  
    self.projects = Table('projects', self.metadata,
        Column('CG_ID_Project', String(15), primary_key=True),
        Column('date_analysis', DateTime),
        Column('organism', String(30)),
        Column('ST', SmallInteger),
      )

    indata = os.listdir(self.config["folders"]["profiles"])
    for file in indata:
      #if not self.engine.dialect.has_table(self.engine, 'profile_{}'.format(file)):
        with open("{}/{}".format(self.config["folders"]["profiles"], file), "r") as fh:
          #Sets profile_* headers
          head = fh.readline()
          head = head.rstrip().split('\t')
          pt = Table('profile_{}'.format(file), self.metadata,
            Column(head[0], SmallInteger, primary_key=True),
            Column(head[1], SmallInteger),
            Column(head[2], SmallInteger),
            Column(head[3], SmallInteger),
            Column(head[3], SmallInteger),
            Column(head[4], SmallInteger),
            Column(head[5], SmallInteger),
            Column(head[6], SmallInteger),
            Column(head[7], SmallInteger),
          )
          self.profiletables[file]=pt

    #Attempts to create tables every time a new object is formed. Maybe not ideal
    self.create_tables()

  def create_tables(self):
      if not self.engine.dialect.has_table(self.engine, 'samples'):
        self.samples.create()
      for k,v in self.profiletables.items():
        if not self.engine.dialect.has_table(self.engine, "profile_{}".format(k)):
          self.profiletables[k].create()
          self.init_profiletable(k, v)         
 
  def add_samplesrecord(self, data_dict):
    #TODO: Dictionary must contain all values found in fields list. Otherwise return error.
    inserter = self.samples.insert()
    #TODO: Checks if entry exists. Give a bucket of errors atm if duplicate record. Should use PK!
    if not self.samples.select((self.samples.c.run == data_dict['run']) & (self.samples.c.contig_name == data_dict['contig_name'])).execute().fetchone():
      inserter.execute(data_dict)
 
  def init_profiletable(self, filename, table):
    data = table.insert()
    linedict = dict.fromkeys(table.c.keys())
    with open("{}/{}".format(self.config["folders"]["profiles"], filename), "r") as fh:
        #Skip header
        head = fh.readline()
        head = head.rstrip().split('\t')
        for line in fh:
            line = line.rstrip().split('\t')
            linedict[head[0]] = line[0]
            linedict[head[1]] = line[1]
            linedict[head[2]] = line[2]
            linedict[head[3]] = line[3]
            linedict[head[4]] = line[4]
            linedict[head[5]] = line[5]
            linedict[head[6]] = line[6]
            linedict[head[7]] = line[7]
            data.execute(linedict)
 
  def get_samplescolumns(self):
    """ Returns a dictionary where each column is a key. Makes re-entry easier """
    return dict.fromkeys(self.samples.c.keys())

  def get_all_samplesrecords(self):
    return self.samples.select().execute().fetchall()

  def st2allele(self, entry):
    """ Takes an entire samples table column and returns the correct allele number """
    #Query correct table
    for k,v in self.profiletables.items():
      if k == entry['organism']:
        return self.session.query(eval("v.c.{}".format(entry['loci']))).filter(v.c.ST==entry['assumed_ST']).scalar()

  def alleles2st(self):
    """Takes an allele list for an organism and calculates most likely ST"""
    pdb.set_trace()
    #Currently can't deal with alternatives
    print("pop")

