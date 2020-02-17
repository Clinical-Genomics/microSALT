#!/usr/bin/env python

import json
import os
import pdb
import pytest
import requests
import time
from unittest.mock import patch

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT import preset_config, logger
from microSALT.cli import root

def create_instance():
  dbm = DB_Manipulator(config=preset_config,log=logger)
  dbm.add_rec({'CG_ID_sample':'AAA1234A1'}, 'Samples')
  dbm.add_rec({'CG_ID_sample':'AAA1234A1', 'loci':'mdh', 'contig_name':'NODE_1'}, 'Seq_types')
  dbm.add_rec({'CG_ID_sample':'AAA1234A1', 'gene':'Type 1', 'instance':'Type 1', 'contig_name':'NODE_1'}, 'Resistances')
  dbm.add_rec({'CG_ID_sample':'AAA1234A1','gene':'Type 1', 'instance':'Type 1', 'contig_name':'NODE_1'}, 'Expacs')

  dbm.add_rec({'CG_ID_project':'AAA1234'}, 'Projects')
  dbm.add_rec({'CG_ID_project':'AAA1234','version':'1'}, 'Reports')

  #dbm.add_rec({'CG_ID_sample':'AAA1234'}, 'Expacs')
  #dbm.add_rec({'CG_ID_sample':'AAA1234'}, 'Collections')
