import collections
import logging
import json
import os
import pathlib
import re
import subprocess
import sys

from flask import Flask
from distutils.sysconfig import get_python_lib

__version__ = '2.9.0'

app = Flask(__name__, template_folder='server/templates')
app.config.setdefault('SQLALCHEMY_DATABASE_URI', 'sqlite:///:memory:')
app.config.setdefault('SQLALCHEMY_BINDS', None)
app.config.setdefault('SQLALCHEMY_TRACK_MODIFICATIONS', False)

#Keep track of microSALT installation
wd=os.path.dirname(os.path.realpath(__file__))

# Load configuration
config = ''
logger = ''
default = os.path.join(os.environ['HOME'], '.microSALT/config.json')

if 'MICROSALT_CONFIG' in os.environ:
  try:
    envvar = os.environ['MICROSALT_CONFIG']
    with open(envvar, 'r') as conf:
      config = json.load(conf)
  except Exception as e:
    print("Config error: {}".format(str(e)))
    pass
elif os.path.exists(default):
  try:
    with open(os.path.abspath(default), 'r') as conf:
      config = json.load(conf)
  except Exception as e:
    print("Config error: {}".format(str(e))) 
    pass

# Config dependent section:
if config != '':
  try:
    #Load flask info
    app.config.update(config['database'])

    #Add extrapaths to config
    config['folders']['expec'] = os.path.abspath(os.path.join(pathlib.Path(__file__).parent.parent, 'unique_references/ExPEC.fsa'))
    #Check if release install exists
    for entry in os.listdir(get_python_lib()):
      if 'microSALT-' in entry:
        config['folders']['expec'] = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'expec/ExPEC.fsa'))
        break
    config['folders']['adapters'] = os.path.abspath(os.path.join(os.path.expandvars('$CONDA_PREFIX'), 'share/trimmomatic-0.39-1/adapters/'))

    #Initialize logger
    logger = logging.getLogger('main_logger')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
    logger.addHandler(ch)

    #Create paths mentioned in config
    for entry in config.keys():
      if entry != '_comment':
        if isinstance(config[entry], str) and '/' in config[entry] and entry not in ['genologics']:
          if not config[entry].startswith('/'):
            sys.exit(-1)
          unmade_fldr = os.path.abspath(config[entry])
          if not pathlib.Path(unmade_fldr).exists():
            os.makedirs(unmade_fldr)
            logger.info("Created path {}".format(unmade_fldr))

        #level two
        elif isinstance(config[entry], collections.Mapping):
          for thing in config[entry].keys():
            if isinstance(config[entry][thing], str) and '/' in config[entry][thing] and entry not in ['genologics']:
              #Special string, mangling
              if thing == 'log_file':
                unmade_fldr = os.path.dirname(config[entry][thing])
                bash_cmd="touch {}".format(config[entry][thing])
                proc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                output, error = proc.communicate()
              elif thing == 'SQLALCHEMY_DATABASE_URI':
                db_file = re.search('sqlite:///(.+)', config['database']['SQLALCHEMY_DATABASE_URI']).group(1)
                unmade_fldr = os.path.dirname(db_file)
                bash_cmd="touch {}".format(db_file) 
                proc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                output, error = proc.communicate()
              else:
                unmade_fldr = config[entry][thing]
              if not pathlib.Path(unmade_fldr).exists():
                os.makedirs(unmade_fldr)
                logger.info("Created path {}".format(unmade_fldr))

    fh = logging.FileHandler(os.path.expanduser(config['folders']['log_file']))
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)

  except Exception as e:
    print("Config error: {}".format(str(e)))
    pass
