# -*- coding: utf-8 -*-

import os
import sys
import yaml
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

__version__ = '1.0.0'

app = Flask(__name__,instance_relative_config=True, template_folder='server/templates')
try:
  app.config.from_pyfile('sqlalchemy_config.py')
except Exception as e:
  print("Unable to load microSALT. Do you have a sqlalchemy_config.py file in folder instance/ ?")
  sys.exit(-1)
db = SQLAlchemy(app)
