#!/usr/bin/env python

import os
import pathlib


# TODO: Add expec loading into the config loading.... however we do that now
def test_existence(exp_config, config):
  """Checks that the configuration contains certain key variables"""
  # level one
  config_level_one = config.keys()
  for entry in exp_config.keys():
    if entry != 'dry':
      assert entry in config_level_one

      # level two
      if isinstance(config[entry], dict):
        config_level_two = config[entry].keys()
        for thing in exp_config[entry]:
          assert thing in config_level_two

def test_reverse_existence(exp_config, config):
  """Check that the configuration doesn't contain outdated variables"""

  # level one
  config_level_one = exp_config.keys()
  for entry in config.keys():
    if entry not in ['_comment']:
      assert entry in config_level_one

      # level two
      config_level_two = exp_config[entry]
      if isinstance(config[entry], dict):
        for thing in config[entry].keys():
          if thing != '_comment':
            assert thing in config_level_two

def test_paths(config):
  """Tests existence for all paths mentioned in variables"""
  # level one
  for entry in config.keys():
    if entry != '_comment':
      if isinstance(config[entry], str) and '/' in config[entry] and entry not in ['database', 'genologics']:
        unmade_fldr = config[entry]
        # Embed logic to expand vars and user here
        unmade_fldr = os.path.expandvars(unmade_fldr)
        unmade_fldr = os.path.expanduser(unmade_fldr)
        unmade_fldr = os.path.abspath(unmade_fldr)
        assert (pathlib.Path(unmade_fldr).exists())
    
      # level two
      elif isinstance(config[entry], dict):
        for thing in config[entry].keys():
          if isinstance(config[entry][thing], str) and '/' in config[entry][thing] and entry not in ['database', 'genologics']:
            unmade_fldr = config[entry][thing]
            # Embed logic to expand vars and user here
            unmade_fldr = os.path.expandvars(unmade_fldr)
            unmade_fldr = os.path.expanduser(unmade_fldr)
            unmade_fldr = os.path.abspath(unmade_fldr)
            assert (pathlib.Path(unmade_fldr).exists())
