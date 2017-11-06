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

  def __init__(self):
    self.engine = create_engine("mysql+pymysql://{}:{}@{}:{}/{}".format(self.mysql['user'],self.mysql['pw'],self.mysql['host'],self.mysql['port'],self.mysql['db']))
    self.metadata = MetaData(self.engine)
    self.blasttable = Table('blasttable', self.metadata,
        Column('run', String(40), primary_key=True),
        Column('date_analysis', DateTime),
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

  def create_blasttable(self):
      if not self.engine.dialect.has_table(self.engine, 'blasttable'):
        self.blasttable.create()
      else:
        print("Blasttable already exists!")
  
  def add_blastrecord(self, data_dict):
    #TODO: Dictionary must contain all values found in fields list. Otherwise return error.
    inserter = self.blasttable.insert()
    #TODO: Checks if entry exists. Give a bucket of errors atm if duplicate record.
    if not self.blasttable.select((self.blasttable.c.run == data_dict['run']) & (self.blasttable.c.loci == data_dict['contig_name'])).execute().fetchone():
      inserter.execute(data_dict)
  
  def get_blastcolumns(self):
    """ Returns a dictionary where each column is a key. Makes re-entry easier """
    return dict.fromkeys(self.blasttable.c.keys())

  def get_all_blastrecords(self):
    return self.blasttable.select().execute().fetchall()

  #def 

