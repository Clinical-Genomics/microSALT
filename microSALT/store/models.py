"""Table definitions for profiles databases. Bit special since it spawns multiple tables.
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import re
from sqlalchemy import *
import yaml

#TODO: Rewrite this implementation into orm format
class Profiles:

  def __init__(self, metadata, config):
    self.tables = dict()

    indata = os.listdir(config["folders"]["profiles"])
    for file in indata:
      #if not self.engine.dialect.has_table(self.engine, 'profile_{}'.format(file)):
        with open("{}/{}".format(config["folders"]["profiles"], file), "r") as fh:
          #Sets profile_* headers
          head = fh.readline()
          head = head.rstrip().split('\t')
          index = 0
          while index < len(head):
            if index == 0:
              header = "Table('profile_{}'.format(file), metadata,".format(file)
              header +="Column(head[{}], SmallInteger, primary_key=True),".format(index)
            else:
              header +="Column(head[{}], SmallInteger),".format(index)
            index = index+1
          header +=")"
          pt = eval(header)
          self.tables[file]=pt

