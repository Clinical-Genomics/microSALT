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

saureus_mlst=[
{'CG_ID_sample':'AAA1234A1',
 'loci':'arcC',
 'allele':'6',
 'contig_name':'NODE_1',
 'contig_length':'592262',
 'contig_coverage':'150',
 'identity':'100.0',
 'evalue':'0.0',
 'bitscore':'900',
 'st_predictor':'1',
 'span':'1.0',
 'subject_length':'456',
 'contig_start':'588316',
 'contig_end':'588771',
},
{'CG_ID_sample':'AAA1234A2',
 'loci':'aroE',
 'allele':'57',
 'contig_name':'NODE_1',
 'contig_length':'592262',
 'contig_coverage':'150',
 'identity':'100.0',
 'evalue':'0.0',
 'bitscore':'900',
 'st_predictor':'1',
 'span':'1.0',
 'subject_length':'456',
 'contig_start':'588316',
 'contig_end':'588771',
},
{'CG_ID_sample':'AAA1234A3',
 'loci':'glpF',
 'allele':'45',
 'contig_name':'NODE_1',
 'contig_length':'592262',
 'contig_coverage':'150',
 'identity':'100.0',
 'evalue':'0.0',
 'bitscore':'900',
 'st_predictor':'1',
 'span':'1.0',
 'subject_length':'456',
 'contig_start':'588316',
 'contig_end':'588771',
},
{'CG_ID_sample':'AAA1234A4',
 'loci':'gmk',
 'allele':'2',
 'contig_name':'NODE_1',
 'contig_length':'592262',
 'contig_coverage':'150',
 'identity':'100.0',
 'evalue':'0.0',
 'bitscore':'900',
 'st_predictor':'1',
 'span':'1.0',
 'subject_length':'456',
 'contig_start':'588316',
 'contig_end':'588771',
},
{'CG_ID_sample':'AAA1234A5',
 'loci':'pta',
 'allele':'7',
 'contig_name':'NODE_1',
 'contig_length':'592262',
 'contig_coverage':'150',
 'identity':'100.0',
 'evalue':'0.0',
 'bitscore':'900',
 'st_predictor':'1',
 'span':'1.0',
 'subject_length':'456',
 'contig_start':'588316',
 'contig_end':'588771',
},
{'CG_ID_sample':'AAA1234A6',
 'loci':'tpi',
 'allele':'58',
 'contig_name':'NODE_1',
 'contig_length':'592262',
 'contig_coverage':'150',
 'identity':'100.0',
 'evalue':'0.0',
 'bitscore':'900',
 'st_predictor':'1',
 'span':'1.0',
 'subject_length':'456',
 'contig_start':'588316',
 'contig_end':'588771',
},
{'CG_ID_sample':'AAA1234A7',
 'loci':'yqiL',
 'allele':'52',
 'contig_name':'NODE_1',
 'contig_length':'592262',
 'contig_coverage':'150',
 'identity':'100.0',
 'evalue':'0.0',
 'bitscore':'900',
 'st_predictor':'1',
 'span':'1.0',
 'subject_length':'456',
 'contig_start':'588316',
 'contig_end':'588771',
}]


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

def add_mlst():
  dbm = DB_Manipulator(config=preset_config,log=logger)
  for entry in saureus_mlst:
    dbm.add_rec(entry, 'Seq_types')
