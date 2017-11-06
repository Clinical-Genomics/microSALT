""" Initial script to deliver and fetch data from a database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pymysql
import re
from sqlalchemy import *
import yaml

import pdb # debug

class DB_Manipulator:

  # Keeping mysql.yml seperate lets us share the main config one without security issues
  with open("{}/mysql.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
    mysql = yaml.load(conf)

  def __init__(self, config):
    self.config = config
    self.engine = create_engine("mysql+pymysql://{}:{}@{}:{}/{}".format(self.mysql['user'],self.mysql['pw'],self.mysql['host'],self.mysql['port'],self.mysql['db']))
    self.metadata = MetaData(self.engine)
    ## Tables
    self.blasttable = Table('blast', self.metadata,
        Column('run', String(40), primary_key=True),
        Column('date_analysis', DateTime),
        Column('organism', String(30)),
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
    #TODO: DAT REDUNDANCY
    self.profiletable = Table('profiles', self.metadata,
        Column('organism', String(30), primary_key=True),
        Column('ST', SmallInteger, primary_key=True),
        Column('loci_1', SmallInteger),
        Column('loci_2', SmallInteger),
        Column('loci_3', SmallInteger),
        Column('loci_4', SmallInteger),
        Column('loci_5', SmallInteger),
        Column('loci_6', SmallInteger),
        Column('loci_7', SmallInteger),
    )
    self.organismloci = Table('organismloci', self.metadata,
        Column('organism', String(30), primary_key=True),
        Column('loci_1', String(5)),
        Column('loci_2', String(5)),
        Column('loci_3', String(5)),
        Column('loci_4', String(5)),
        Column('loci_5', String(5)),
        Column('loci_6', String(5)),
        Column('loci_7', String(5)),
    )

    #Attempts to create tables every time a new object is formed. Maybe not ideal
    self.create_tables()

  def create_tables(self):
      if not self.engine.dialect.has_table(self.engine, 'blast'):
        self.blasttable.create()
      if not self.engine.dialect.has_table(self.engine, 'organismloci') and not self.engine.dialect.has_table(self.engine, 'profiles'):
        self.profiletable.create()
        self.organismloci.create()
        self.init_profilerecords()

  def add_blastrecord(self, data_dict):
    #TODO: Dictionary must contain all values found in fields list. Otherwise return error.
    inserter = self.blasttable.insert()
    #TODO: Checks if entry exists. Give a bucket of errors atm if duplicate record. Should use PK!
    if not self.blasttable.select((self.blasttable.c.run == data_dict['run']) & (self.blasttable.c.loci == data_dict['contig_name'])).execute().fetchone():
      inserter.execute(data_dict)
 
  def init_profilerecords(self):
    indata = os.listdir(self.config["folders"]["profiles"])
    data = self.profiletable.insert()
    headers = self.organismloci.insert()

    headdict = self.get_organismlocicolumns()
    linedict = self.get_profilecolumns()

    for file in indata:
      headdict['organism'] = file
      linedict['organism'] = file
      with open("{}/{}".format(self.config["folders"]["profiles"], file), "r") as fh:
        #Send first row to organismloci table
        head = fh.readline()
        head = head.rstrip().split('\t')
        headdict['loci_1'] = head[1]
        headdict['loci_2'] = head[2]
        headdict['loci_3'] = head[3]
        headdict['loci_4'] = head[4]
        headdict['loci_5'] = head[5]
        headdict['loci_6'] = head[6]
        headdict['loci_7'] = head[7]
        headers.execute(headdict)

        #Rest to profiles table
        for line in fh:
          line = line.rstrip().split('\t')
          linedict['ST'] = line[0]
          linedict['loci_1'] = line[1]
          linedict['loci_2'] = line[2]
          linedict['loci_3'] = line[3]
          linedict['loci_4'] = line[4]
          linedict['loci_5'] = line[5]
          linedict['loci_6'] = line[6]
          linedict['loci_7'] = line[7]
          data.execute(linedict)
 
  def get_blastcolumns(self):
    """ Returns a dictionary where each column is a key. Makes re-entry easier """
    return dict.fromkeys(self.blasttable.c.keys())

  def get_profilecolumns(self):
    return dict.fromkeys(self.profiletable.c.keys())

  def get_organismlocicolumns(self):
    return dict.fromkeys(self.organismloci.c.keys())

  def get_all_blastrecords(self):
    return self.blasttable.select().execute().fetchall()

