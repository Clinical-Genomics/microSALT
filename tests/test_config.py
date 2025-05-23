#!/usr/bin/env python

import collections
import os
import pathlib
import pytest

from microSALT import preset_config

@pytest.fixture
def exp_config():
  precon = {
    'slurm_header': {'time', 'threads', 'qos', 'job_prefix', 'project', 'type'},
    'regex': {'file_pattern', 'mail_recipient', 'verified_organisms'},
    'folders': {'results', 'reports', 'log_file', 'seqdata', 'profiles', 'references', 'resistances', 'genomes', 'expec', 'adapters', 'credentials'},
    'threshold': {'mlst_id', 'mlst_novel_id', 'mlst_span', 'motif_id', 'motif_span', 'total_reads_warn', 'total_reads_fail',
                  'NTC_total_reads_warn', 'NTC_total_reads_fail', 'mapped_rate_warn', 'mapped_rate_fail', 'duplication_rate_warn',
                  'duplication_rate_fail', 'insert_size_warn', 'insert_size_fail', 'average_coverage_warn', 'average_coverage_fail',
                  'bp_10x_warn', 'bp_10x_fail', 'bp_30x_warn', 'bp_50x_warn', 'bp_100x_warn'},
    'database': {'SQLALCHEMY_DATABASE_URI', 'SQLALCHEMY_TRACK_MODIFICATIONS', 'DEBUG'},
    'genologics': {'baseuri', 'username', 'password'},
    'pubmlst': {'client_id', 'client_secret'},
    'pasteur': {'client_id', 'client_secret'},
    'dry': True,
  }
  return precon

def test_existence(exp_config):
  """Checks that the configuration contains certain key variables"""
  # level one
  config_level_one = preset_config.keys()
  for entry in exp_config.keys():
    if entry != 'dry':
      assert entry in config_level_one

      # level two
      if isinstance(preset_config[entry], collections.Mapping):
        config_level_two = preset_config[entry].keys()
        for thing in exp_config[entry]:
          assert thing in config_level_two

def test_reverse_existence(exp_config):
  """Check that the configuration doesn't contain outdated variables"""

  # level one
  config_level_one = exp_config.keys()
  for entry in preset_config.keys():
    if entry not in ['_comment']:
      assert entry in config_level_one

      # level two
      config_level_two = exp_config[entry]
      if isinstance(preset_config[entry], collections.Mapping):
        for thing in preset_config[entry].keys():
          if thing != '_comment':
            assert thing in config_level_two

def test_paths(exp_config):
  """Tests existence for all paths mentioned in variables"""
  # level one
  for entry in preset_config.keys():
    if entry != '_comment':
      if isinstance(preset_config[entry], str) and '/' in preset_config[entry] and entry not in ['database', 'genologics']:
        unmade_fldr = preset_config[entry]
        # Embed logic to expand vars and user here
        unmade_fldr = os.path.expandvars(unmade_fldr)
        unmade_fldr = os.path.expanduser(unmade_fldr)
        unmade_fldr = os.path.abspath(unmade_fldr)
        assert (pathlib.Path(unmade_fldr).exists())
    
      # level two
      elif isinstance(preset_config[entry], collections.Mapping):
        for thing in preset_config[entry].keys():
          if isinstance(preset_config[entry][thing], str) and '/' in preset_config[entry][thing] and entry not in ['database', 'genologics']:
            unmade_fldr = preset_config[entry][thing]
            # Embed logic to expand vars and user here
            unmade_fldr = os.path.expandvars(unmade_fldr)
            unmade_fldr = os.path.expanduser(unmade_fldr)
            unmade_fldr = os.path.abspath(unmade_fldr)
            assert (pathlib.Path(unmade_fldr).exists())
