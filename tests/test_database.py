#!/usr/bin/env python

import json
import os
import pathlib
import pdb
import pytest
import re
import requests
import sys
import time

from distutils.sysconfig import get_python_lib
from unittest.mock import patch

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT import preset_config, logger
from microSALT.cli import root

@pytest.fixture
def mlst_data():
  pfile = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_mlst.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      pfile = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_mlst.json'))
  try:
    with open(pfile) as json_file:
      data = json.load(json_file)
  except Exception as e:
    print("Unable to read provided sample info file as json. Exiting..")
    sys.exit(-1)
  return data

@pytest.fixture
def dbm():
  db_file = re.search('sqlite:///(.+)', preset_config['database']['SQLALCHEMY_DATABASE_URI']).group(1)
  os.remove(db_file)
  dbm = DB_Manipulator(config=preset_config,log=logger)
  dbm.create_tables()
  return dbm


def test_create_every_table(dbm):
  assert dbm.engine.dialect.has_table(dbm.engine, 'samples')
  assert dbm.engine.dialect.has_table(dbm.engine, 'seq_types')
  assert dbm.engine.dialect.has_table(dbm.engine, 'resistances')
  assert dbm.engine.dialect.has_table(dbm.engine, 'expacs')
  assert dbm.engine.dialect.has_table(dbm.engine, 'projects')
  assert dbm.engine.dialect.has_table(dbm.engine, 'reports')
  assert dbm.engine.dialect.has_table(dbm.engine, 'collections')


def test_add_rec(dbm):
  #Adds records to all databases
  dbm.add_rec({'ST':'130','arcC':'6','aroE':'57','glpF':'45','gmk':'2','pta':'7','tpi':'58','yqiL':'52','clonal_complex':'CC1'}, dbm.profiles['staphylococcus_aureus'])
  assert len(dbm.query_rec(dbm.profiles['staphylococcus_aureus'], {'ST':'130'})) == 1
  assert len(dbm.query_rec(dbm.profiles['staphylococcus_aureus'], {'ST':'-1'})) == 0
  dbm.add_rec({'ST':'130','arcC':'6','aroE':'57','glpF':'45','gmk':'2','pta':'7','tpi':'58','yqiL':'52','clonal_complex':'CC1'}, dbm.novel['staphylococcus_aureus'])
  assert len(dbm.query_rec(dbm.novel['staphylococcus_aureus'], {'ST':'130'})) == 1
  assert len(dbm.query_rec(dbm.novel['staphylococcus_aureus'], {'ST':'-1'})) == 0


  dbm.add_rec({'CG_ID_sample':'ADD1234A1'}, 'Samples')
  assert len(dbm.query_rec('Samples', {'CG_ID_sample':'ADD1234A1'})) > 0
  assert len(dbm.query_rec('Samples', {'CG_ID_sample':'XXX1234A10'})) == 0

  dbm.add_rec({'CG_ID_sample':'ADD1234A1', 'loci':'mdh', 'contig_name':'NODE_1'}, 'Seq_types')
  assert len(dbm.query_rec('Seq_types', {'CG_ID_sample':'ADD1234A1', 'loci':'mdh', 'contig_name':'NODE_1'})) > 0
  assert len(dbm.query_rec('Seq_types', {'CG_ID_sample':'XXX1234A10', 'loci':'mdh', 'contig_name':'NODE_1'})) == 0

  dbm.add_rec({'CG_ID_sample':'ADD1234A1', 'gene':'Type 1', 'instance':'Type 1', 'contig_name':'NODE_1'}, 'Resistances')
  assert len(dbm.query_rec('Resistances',{'CG_ID_sample':'ADD1234A1', 'gene':'Type 1', 'instance':'Type 1', 'contig_name':'NODE_1'})) > 0
  assert len(dbm.query_rec('Resistances',{'CG_ID_sample':'XXX1234A10', 'gene':'Type 1', 'instance':'Type 1', 'contig_name':'NODE_1'})) == 0

  dbm.add_rec({'CG_ID_sample':'ADD1234A1','gene':'Type 1', 'instance':'Type 1', 'contig_name':'NODE_1'}, 'Expacs')
  assert len(dbm.query_rec('Expacs',{'CG_ID_sample':'ADD1234A1','gene':'Type 1', 'instance':'Type 1', 'contig_name':'NODE_1'})) > 0
  assert len(dbm.query_rec('Expacs',{'CG_ID_sample':'XXX1234A10','gene':'Type 1', 'instance':'Type 1', 'contig_name':'NODE_1'})) == 0

  dbm.add_rec({'CG_ID_project':'ADD1234'}, 'Projects')
  assert len(dbm.query_rec('Projects',{'CG_ID_project':'ADD1234'})) > 0 
  assert len(dbm.query_rec('Projects',{'CG_ID_project':'XXX1234'})) == 0

  dbm.add_rec({'CG_ID_project':'ADD1234','version':'1'}, 'Reports')
  assert len(dbm.query_rec('Reports',{'CG_ID_project':'ADD1234','version':'1'})) > 0
  assert len(dbm.query_rec('Reports',{'CG_ID_project':'XXX1234','version':'1'})) == 0

  dbm.add_rec({'CG_ID_sample':'ADD1234', 'ID_collection':'MyCollectionFolder'}, 'Collections')
  assert len(dbm.query_rec('Collections',{'CG_ID_sample':'ADD1234', 'ID_collection':'MyCollectionFolder'})) > 0 
  assert len(dbm.query_rec('Collections',{'CG_ID_sample':'XXX1234', 'ID_collection':'MyCollectionFolder'})) == 0

def test_upd_rec(dbm):
  dbm.add_rec({'CG_ID_sample':'UPD1234A1'}, 'Samples')
  assert len(dbm.query_rec('Samples', {'CG_ID_sample':'UPD1234A1'})) == 1
  assert len(dbm.query_rec('Samples', {'CG_ID_sample':'UPD1234A2'})) == 0
  dbm.upd_rec({'CG_ID_sample':'UPD1234A1'}, 'Samples', {'CG_ID_sample':'UPD1234A2'})
  assert len(dbm.query_rec('Samples', {'CG_ID_sample':'UPD1234A1'})) == 0
  assert len(dbm.query_rec('Samples', {'CG_ID_sample':'UPD1234A2'})) == 1

def test_allele_ranker(dbm, mlst_data):
  for entry in mlst_data:
    dbm.add_rec(entry, 'Seq_types')
  dbm.add_rec({'CG_ID_sample':'MLS1234A1', 'CG_ID_project':'MLS1234','organism':'staphylococcus_aureus'}, 'Samples')
  assert dbm.alleles2st('MLS1234A1') == 130
  best_alleles = {'arcC': {'contig_name': 'NODE_1', 'allele': 6}, 'aroE': {'contig_name': 'NODE_1', 'allele': 57}, 'glpF': {'contig_name': 'NODE_1', 'allele': 45}, 'gmk': {'contig_name': 'NODE_1', 'allele': 2}, 'pta': {'contig_name': 'NODE_1', 'allele': 7}, 'tpi': {'contig_name': 'NODE_1', 'allele': 58}, 'yqiL': {'contig_name': 'NODE_1', 'allele': 52}}
  assert dbm.bestAlleles('MLS1234A1') == best_alleles

  for entry in mlst_data:
    entry['allele'] = 0
    entry['CG_ID_sample'] = 'MLS1234A2'
    dbm.add_rec(entry, 'Seq_types') 
  dbm.alleles2st('MLS1234A2') == -1

def test_get_and_set_report(dbm):
  dbm.add_rec({'CG_ID_sample':'ADD1234A1', 'method_sequencing':'1000:1'}, 'Samples')
  dbm.add_rec({'CG_ID_project':'ADD1234','version':'1'}, 'Reports')
  assert dbm.get_report('ADD1234').version == 1
  dbm.upd_rec({'CG_ID_sample':'ADD1234A1', 'method_sequencing':'1000:1'}, 'Samples', {'CG_ID_sample':'ADD1234A1', 'method_sequencing':'1000:2'})
  dbm.set_report('ADD1234')
  assert dbm.get_report('ADD1234').version != 1
