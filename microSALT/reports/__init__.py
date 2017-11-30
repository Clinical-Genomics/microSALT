import os
import yaml
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

__version__ = '1.0.0'

app = Flask(__name__,instance_relative_config=True)
app.config.from_pyfile('config.py')
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)
