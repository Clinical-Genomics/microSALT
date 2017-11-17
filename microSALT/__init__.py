import os
import yaml
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

__version__ = '1.0.0'

with open("{}/configs/mysql.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
  mysql = yaml.load(conf)

app = Flask(__name__)
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SQLALCHEMY_DATABASE_URI'] = "mysql+pymysql://{}:{}@{}:{}/{}".format(mysql['user'],mysql['pw'],mysql['host'],mysql['port'],mysql['db'])

db = SQLAlchemy(app)


