#!/usr/bin/env python

import datetime
import glob
import json
import logging
import os
import pathlib
import pdb
import pytest
import re
import sys

from distutils.sysconfig import get_python_lib

from microSALT import preset_config, logger
from microSALT.utils.reporter import Reporter
from microSALT.utils.referencer import Referencer
from microSALT.store.db_manipulator import DB_Manipulator

def date_hook(json_dict):
    for (key, value) in json_dict.items():
        try:
            json_dict[key] = datetime.datetime.strptime(value, "%Y-%m-%dT%H:%M:%S")
        except:
            pass
    return json_dict

@pytest.fixture
def reference_data():
  testdata = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_samples.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      testdata = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_samples.json'))
  with open(testdata) as json_file:
    data = json.load(json_file)
  return data

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
def resistance_data():
  pfile = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_resistance.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()): 
    if 'microSALT-' in entry:
      pfile = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_resistance.json'))
  try:
    with open(pfile) as json_file:
      data = json.load(json_file)
  except Exception as e:
    print("Unable to read provided sample info file as json. Exiting..")
    sys.exit(-1)
  return data

@pytest.fixture
def expec_data():
  pfile = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_expec.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()): 
    if 'microSALT-' in entry:
      pfile = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_expec.json'))
  try:
    with open(pfile) as json_file:
      data = json.load(json_file)
  except Exception as e:
    print("Unable to read provided sample info file as json. Exiting..")
    sys.exit(-1)
  return data

@pytest.fixture
def report_data():
  pfile = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_reports.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      pfile = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_reports.json'))
  try:
    with open(pfile) as json_file:
      data = json.load(json_file)
      #xata = json.load(json_file)
      #data = list()
      #for entry in xata:
      #  data.append( json.loads(json.dumps(entry), object_hook=date_hook ))
  except Exception as e:
    print("Unable to read provided sample info file as json. Exiting..")
    sys.exit(-1)
  return data

@pytest.fixture
def project_data():
  pfile = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'tests/testdata/sampleinfo_projects.json'))
  #Check if release install exists
  for entry in os.listdir(get_python_lib()):
    if 'microSALT-' in entry:
      pfile = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'testdata/sampleinfo_projects.json'))
  try:
    with open(pfile) as json_file:
      data = json.load(json_file, object_hook=date_hook)
  except Exception as e:
    print("Unable to read provided sample info file as json. Exiting..")
    sys.exit(-1)
  return data

@pytest.fixture
def dbm():
  db_file = re.search('sqlite:///(.+)', preset_config['database']['SQLALCHEMY_DATABASE_URI']).group(1)
  dbm = DB_Manipulator(config=preset_config,log=logger)
  dbm.create_tables()

  assert dbm.engine.dialect.has_table(dbm.engine, 'samples')
  assert dbm.engine.dialect.has_table(dbm.engine, 'seq_types')
  assert dbm.engine.dialect.has_table(dbm.engine, 'resistances')
  assert dbm.engine.dialect.has_table(dbm.engine, 'expacs')
  assert dbm.engine.dialect.has_table(dbm.engine, 'projects')
  assert dbm.engine.dialect.has_table(dbm.engine, 'reports')
  assert dbm.engine.dialect.has_table(dbm.engine, 'collections')
  return dbm

@pytest.fixture
def init_database(dbm, mlst_data, resistance_data, expec_data, report_data, project_data):
  for xentry in project_data:
    dbm.add_rec(xentry, 'Projects')
  for entry in mlst_data:
    dbm.add_rec(entry, 'Seq_types')
  for bentry in resistance_data:
    dbm.add_rec(bentry, 'Resistances')
  for centry in expec_data:
    dbm.add_rec(centry, 'Expacs')
  for dentry in report_data:
    dbm.add_rec(dentry, 'Reports')

@pytest.fixture
def init_references(reference_data):
  ref_obj = Referencer(config=preset_config, log=logger, sampleinfo=reference_data)
  ref_obj.identify_new(reference_data[0].get('CG_ID_project'),project=True)
  ref_obj.update_refs()

@pytest.fixture
def reporter(reference_data):
  reporter_obj = Reporter(config=preset_config, log=logger, sampleinfo=reference_data[0], name="MIC1234A1", output="/tmp/MLST")
  return reporter_obj

def test_motif(init_database, reporter):
  reporter.create_subfolders()
  reporter.gen_motif(motif="resistance")
  assert len( glob.glob("{}/AAA1234_resistance*".format(reporter.output))) > 0

  reporter.gen_motif(motif="expec")
  assert len( glob.glob("{}/AAA1234_expec*".format(reporter.output))) > 0

def test_deliveryreport(init_database, reporter):
  reporter.create_subfolders()
  reporter.gen_delivery()
  assert len( glob.glob("{}/deliverables/999999_deliverables.yaml".format(preset_config['folders']['reports']))) > 0

def test_jsonreport(init_database, reporter):
  reporter.create_subfolders()
  reporter.gen_json()
  assert len( glob.glob("{}/json/AAA1234.json".format(preset_config['folders']['reports']))) > 0
