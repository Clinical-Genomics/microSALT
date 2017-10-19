"""This initial script creates sbatch jobs
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import re
import time
import yaml

class Job_Creator():

  fileformat = re.compile('(\d{1}_\d{6}_\w{9}_.{10,12}_\w{8,12}_)(\d{1})(.fastq.gz)') 
  with open("{}/config.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
      config = yaml.load(conf)

  def __init__(self, indir):
    self.now = time.strftime("%Y.%m.%d_%H.%M.%S")
    self.indir = os.path.abspath(indir)
    self.name = os.path.basename(os.path.normpath(indir))
    self.batchfile = ""
    self.outdir = ""

  def create_header(self):
    batchfile = open(self.batchfile, "w+")
    batchfile.write("#!/bin/bash -l\n\n")
    batchfile.write("#SBATCH -A {}\n".format(self.config["slurm_header"]["project"]))
    batchfile.write("#SBATCH -p {}\n".format(self.config["slurm_header"]["type"]))
    batchfile.write("#SBATCH -n {}\n".format(self.config["slurm_header"]["threads"]))
    batchfile.write("#SBATCH -t {}\n".format(self.config["slurm_header"]["time"]))
    batchfile.write("#SBATCH -J {}_{}\n".format(self.config["slurm_header"]["job_name"], self.now))
    batchfile.write("#SBATCH --qos {}\n\n".format(self.config["slurm_header"]["qos"]))
    batchfile.close()

  def create_spadesjob(self):
    ## Write intial cli
    batchfile = open(self.batchfile, "a+")
    #memory is actually 128 per node regardless of cores.
    batchfile.write("spades.py --threads {} --memory {} -o {}/assembly"\
    .format(self.config["slurm_header"]["threads"], 8*int(self.config["slurm_header"]["threads"]), self.outdir))

    ## Establish valid file pairs
    files = os.listdir(self.indir)
    verified_files = list()
    while len(files) > 0:
      file_parts = self.fileformat.match( files.pop(0) )
      #If file meets standard format, find pair
      if file_parts:
        if file_parts[2] == '1':
          pairno = '2'
        elif file_parts[2] == '2':
          pairno = '1'
        else:
          print("No pair mate found. Add an error message.")

        pairname = "{}{}{}".format(file_parts[1], pairno, file_parts[3])
        if pairname in files:
          files.pop( files.index(pairname) )
          verified_files.append(file_parts[0])
          verified_files.append(pairname)
        else:
          print("No pair mate found. Add an error message.")

    ## Write valid file pairs
    pairno = 0
    for file in verified_files:
      #If index is even, raise pair no
      if (verified_files.index(file) & 1) == 0:
        pairno = pairno + 1
      #Set read no
      readno = verified_files.index(file) % 2 + 1
      batchfile.write(" --pe{}-{} {}/{}".format(pairno, readno, self.indir, file))
    batchfile.write("\n\n")
    batchfile.close()

  def create_blastjob(self):
    #index database
    batchfile = open(self.batchfile, "a+")
    #TODO: Organism needs to be established in the future
    batchfile.write("cd {} && makeblastdb -in {}/{}.xmfa -dbtype nucl -parse_seqids -out {}\n".format(\
    self.config["folders"]["references"], self.config["folders"]["references"], self.config["organism"], self.config["organism"]))
    #create run
    blast_format = "7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send"
    batchfile.write("blastn -db {}/{} -query {}/assembly/contigs.fasta -out {}/loci_query_tab.txt -num_threads {} -max_target_seqs 1 -outfmt {}\n\n".format(\
    self.config["folders"]["references"], self.config["organism"], self.outdir, self.outdir, self.config["slurm_header"]["threads"], blast_format))


  def create_job(self):
    self.outdir = "{}/{}_{}".format(self.config["folders"]["results"],self.name, self.now)
    if not os.path.exists(self.outdir):
      os.makedirs(self.outdir)
    self.batchfile = "{}/runfile.sbatch".format(self.outdir)

    self.create_header()
    self.create_spadesjob()
    self.create_blastjob()
