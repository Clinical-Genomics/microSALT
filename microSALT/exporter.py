import click
import os
import re
import time
import yaml

from microSALT import db_manipulator

import pdb # debug

class Exporter():

  def __init__(self, config):
    self.config = config
    self.db_manip=db_manipulator.DB_Manipulator()

  def std_export(self):
    now = time.strftime("%Y.%m.%d_%H.%M.%S")
    header = self.db_manip.get_blastcolumns()
    data = self.db_manip.get_all_blastrecords()
    data.insert(0, header)
    outfile = open("blast_{}.dump".format(now), "w+")
    for row in data:
      line = "|"
      for field in row:
        line = line + "{}|".format(str(field))
      outfile.write("{}\n".format(line))
