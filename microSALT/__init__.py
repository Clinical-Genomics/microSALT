import os
import yaml
from flask import Flask
from flask_sqlalchemy import SQLAlchemy

with open("{}/configs/mysql.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
  mysql = yaml.load(conf)

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = "mysql+pymysql://{}:{}@{}:{}/{}".format(mysql['user'],mysql['pw'],mysql['host'],mysql['port'],mysql['db'])
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)


__version__ = '1.0.0'
