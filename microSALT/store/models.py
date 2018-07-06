"""Table definitions for profiles databases. Bit special since it spawns multiple tables.
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import os
from sqlalchemy import *

#TODO: Rewrite this implementation into orm format
class Profiles:

  def __init__(self, metadata, config):
    self.tables = dict()
    self.metadata = metadata
    self.config = config
    indata = os.listdir(self.config["folders"]["profiles"])
    for file in indata:
      self.add_table(file)

  def add_table(self, file):
    with open("{}/{}".format(self.config["folders"]["profiles"], file), "r") as fh:
      #Sets profile_* headers
      head = fh.readline()
      head = head.rstrip().split('\t')
      index = 0
      header = "Table('profile_{}'.format(file), self.metadata,".format(file)
      while index < len(head):
        # Set ST as PK
        if index == 0:
          header +="Column(head[{}], SmallInteger, primary_key=True),".format(index)
        # Set Clonal complex as string
        elif index == len(head)-1:
          header +="Column(head[{}], String(16)),".format(index)
        else:
          header +="Column(head[{}], SmallInteger),".format(index)
        index = index+1
      header +=")" 
      pt = eval(header)
      self.tables[file]=pt

