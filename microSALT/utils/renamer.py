#!/usr/bin/env python

import os
import re
from microSALT.store.lims_fetcher import LIMS_Fetcher

class Renamer():

  def __init__(self, infolder, config, log):
    self.config = config
    self.logger = log
    self.infolder = os.path.abspath(infolder)
    self.sampledir = ""
    self.limsfetcher = LIMS_Fetcher(self.logger, self.config)

  def rename_project(self):
    projname = os.path.basename(os.path.normpath(self.infolder))
    fileformat = re.compile(self.config['regex']['rename_pattern'])
  
    for root, dir, files in os.walk(self.infolder):
      dir_parts = os.path.split(root)
      if dir_parts[1] == projname:
        continue
      else:
        samples = self.limsfetcher.lims.get_samples(name=dir_parts[1])
        if samples == []:
          self.logger.warn("Unable to find samples with external name '{}' of {}. Input files already contain internal ids".format(dir_parts[1], projname))
        else:
          for sample in samples:
            if sample.project.id == projname:
              break
          #Change folder
          newfolder = "{}/{}".format(dir_parts[0], sample.id)
          os.rename(root, newfolder)
          #Change files
          for f in files:
            fm = fileformat.match(f)
            id_match = re.search(fm[1], f)
            os.rename("{}/{}".format(newfolder, f), "{}/{}{}{}".format(newfolder, f[:id_match.start()] , sample.id, f[id_match():]) )
        
