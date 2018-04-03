import json
import os
import sys

from flask import Flask

__version__ = '2.0.0'
app = Flask(__name__, template_folder='server/templates')
app.config.setdefault('SQLALCHEMY_DATABASE_URI', 'sqlite:///:memory:')
app.config.setdefault('SQLALCHEMY_BINDS', None)
app.config.setdefault('SQLALCHEMY_TRACK_MODIFICATIONS', False)

# Load configuration
config = ''
defaulto = os.path.join(os.environ['HOME'], '.microSALT/config.json')
if os.path.exists(defaulto):
  try:
    with open(defaulto, 'r') as conf:
      config = json.load(conf)
  except Exception as e:
    pass
elif 'MICROSALT_CONFIG' in os.environ:
  try:
    envvar = os.environ['MICROSALT_CONFIG']
    with open(envvar, 'r') as conf:
      config = json.load(conf)
  except Exception as e:
    pass
# Load flask instance
if config != '':
  try:
    app.config.update(config['database'])
  except Exception as e:
    pass
