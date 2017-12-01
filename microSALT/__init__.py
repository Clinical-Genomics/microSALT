# -*- coding: utf-8 -*-

import os
import yaml
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

__version__ = '1.0.0'

app = Flask(__name__,instance_relative_config=True, template_folder='server/templates')
app.config.from_pyfile('config.py')
db = SQLAlchemy(app)
