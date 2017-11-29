"""This initial script runs jobs
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pdb
import re
import subprocess
import sys
import time
import yaml

from microSALT.job_creator import Job_Creator

class Manager():

  def __init__(self, indir, config, log):
    self.config = config
    self.now = time.strftime("%Y.%m.%d_%H.%M.%S")
    self.indir = os.path.abspath(indir)
    self.name = os.path.basename(os.path.normpath(indir))
    self.logger = log
    self.concat_batch = "{}/{}.sbatch".format(os.getcwd(), self.name)

  def start_analysis(self):
    batchfile = open(self.batchfile, 'w+') 
    #Write an sbatch job per sample
    for (dirpath, dirnames, filenames) in os.walk(self.indir):
      for dir in dirnames:
        organism = subprocess.Popen("cglims get {} Strain".format(dir).split(), stdout=subprocess.PIPE).communicate()
        new_job = Job_Creator(str(organism), self.config, self.log)
        new_job.find_reference_loose()
        new_job.create_job()
        outfile = new_job.get_sbatch()
        #In the future, rather than running the files, just start everything in this list.
        #Make this a dry-run command instead.
        batchfile.write("sbatch {}\n".format(outfile))

    batchfile.close()
