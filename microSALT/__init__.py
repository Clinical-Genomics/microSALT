import json
import os
import sys

from flask import Flask

__version__ = '2.8.15'

app = Flask(__name__, template_folder='server/templates')
app.config.setdefault('SQLALCHEMY_DATABASE_URI', 'sqlite:///:memory:')
app.config.setdefault('SQLALCHEMY_BINDS', None)
app.config.setdefault('SQLALCHEMY_TRACK_MODIFICATIONS', False)

#Keep track of microSALT installation
wd=os.path.dirname(os.path.realpath(__file__))

# Load configuration
config = ''
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
# Load flask instance
if config != '':
  try:
    app.config.update(config['database'])
  except Exception as e:
    print("Config error: {}".format(str(e)))
    pass

