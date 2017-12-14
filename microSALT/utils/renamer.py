#!/usr/bin/env python

import os
import re
from genologics.lims import Lims
# Should probably call these items directly since we're now up to 3 config files
from genologics.config import BASEURI,USERNAME,PASSWORD
from genologics.entities import Project, Sample

#TODO: Add functions to lims_fetcher class. Put fileformat under headers config file.
class Renamer():

  def __init__(self, infolder, config, log):
    self.config = config
    self.logger = log
    self.infolder = os.path.abspath(infolder)
    self.sampledir = ""

  def rename_project(self):
    self.lims = Lims(BASEURI, USERNAME, PASSWORD)
    fileformat = re.compile('(\d{1}_\d{6}_\w{9}_|\d{6}_\w{9}_).{3,12}(_\w{8,12}_\d{1}.fastq.gz|_\d{1}.fastq.gz)')
    projname = os.path.basename(os.path.normpath(self.infolder))
  
    for root, dir, files in os.walk(self.infolder):
      dir_parts = os.path.split(root)
      if dir_parts[1] == projname:
        continue
      else:
        samples = self.lims.get_samples(name=dir_parts[1])
        for sample in samples:
          if sample.project.id == projname:
            break
        #Change folder
        newfolder = "{}/{}".format(dir_parts[0], sample.id)
        os.rename(root, newfolder)
        #Change files
        for f in files:
          filematch = fileformat.match(f)
          os.rename("{}/{}".format(newfolder, f), "{}/{}{}{}".format(newfolder, filematch[1], sample.id, filematch[2]))
