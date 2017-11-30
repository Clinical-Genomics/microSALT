# -*- coding: utf-8 -*-

import os
import yaml
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

__version__ = '1.0.0'

app = Flask(__name__,instance_relative_config=True, template_folder='reports/templates')
app.config.from_pyfile('config.py')
db = SQLAlchemy(app)

## Legacy implementation. Kept just to make the config file more obvious
#with open("{}/config/mysql.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
#  mysql = yaml.load(conf)
#app = Flask(__name__, template_folder='reports/templates')
#app.config['SQLALCHEMY_DATABASE_URI'] = "mysql+pymysql://{}:{}@{}:{}/{}".format(mysql['user'],mysql['pw'],mysql['host'],mysql['port'],mysql['db'])
#app.config['SQLALCHEMY_TRACK_MODIFICATIONS']= False
#app.config['DEBUG']= True
