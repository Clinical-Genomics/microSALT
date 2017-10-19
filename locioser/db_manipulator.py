""" Initial script to deliver and fetch data from a database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pymysql
import re
from sqlalchemy import *
import yaml

import pdb # debug

with open("{}/mysql.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
  mysql = yaml.load(conf)

pdb.set_trace()

engine = create_engine("mysql+pymysql://{}:{}@{}:{}/{}".format(mysql['user'],mysql['pw'],mysql['host'],mysql['port'],mysql['db']))
engine.echo = True
metadata = MetaData(engine)

table1 = Table('testtable', metadata,
  Column('testvalue', Integer),
)
table1.create()

inserter = table1.insert()
inserter.execute(testvalue=9)

s = table1.select()
rs = s.execute()
row = rs.fetchone()
