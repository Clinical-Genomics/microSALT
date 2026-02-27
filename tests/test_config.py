#!/usr/bin/env python

import collections
import collections.abc
import os
import pathlib

from microSALT import preset_config


def test_existence(exp_config):
  """Checks that the configuration contains certain key variables"""
  # level one
  config_level_one = preset_config.keys()
  for entry in exp_config.keys():
    assert entry in config_level_one

    # level two
    if isinstance(preset_config[entry], collections.abc.Mapping):
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
      if isinstance(preset_config[entry], collections.abc.Mapping):
        for thing in preset_config[entry].keys():
          if thing != '_comment':
            assert thing in config_level_two

def test_paths(exp_config):
  """Tests existence for all paths mentioned in variables"""
  # level one
  for entry in preset_config.keys():
    if entry != '_comment':
      if isinstance(preset_config[entry], str) and '/' in preset_config[entry] and entry not in ['database']:
        unmade_fldr = preset_config[entry]
        # Embed logic to expand vars and user here
        unmade_fldr = os.path.expandvars(unmade_fldr)
        unmade_fldr = os.path.expanduser(unmade_fldr)
        unmade_fldr = os.path.abspath(unmade_fldr)
        assert (pathlib.Path(unmade_fldr).exists())
    
      # level two
      elif isinstance(preset_config[entry], collections.abc.Mapping):
        for thing in preset_config[entry].keys():
          if isinstance(preset_config[entry][thing], str) and '/' in preset_config[entry][thing] and entry not in ['database']:
            unmade_fldr = preset_config[entry][thing]
            # Embed logic to expand vars and user here
            unmade_fldr = os.path.expandvars(unmade_fldr)
            unmade_fldr = os.path.expanduser(unmade_fldr)
            unmade_fldr = os.path.abspath(unmade_fldr)
            assert (pathlib.Path(unmade_fldr).exists())
