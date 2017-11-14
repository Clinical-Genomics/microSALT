"""This initial script creates sbatch jobs
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pdb
import re
import sys
import time
import yaml

class Job_Creator():

  fileformat = re.compile('(\d{1}_\d{6}_\w{9}_.{5,12}_\w{8,12}_)(\d{1})(.fastq.gz)') 

  def __init__(self, indir, organism, config, log):
    self.config = config
    self.now = time.strftime("%Y.%m.%d_%H.%M.%S")
    self.indir = os.path.abspath(indir)
    self.name = os.path.basename(os.path.normpath(indir))
    self.organism = organism
    self.logger = log
    self.batchfile = ""
    self.outdir = "{}/{}_{}".format(self.config["folders"]["results"],self.name, self.now)
    self.trimmed_files = dict()

  def verify_fastq(self):
    """ Uses arg indir to return a list of PE fastq tuples fulfilling naming convention """
    files = os.listdir(self.indir)
    if files == []:
      self.logger.error("No fastq files found in specified directory. Exited.")
      sys.exit()
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
          self.logger.error("Some fastq files in directory have no mate. Exited.")
          sys.exit()
        pairname = "{}{}{}".format(file_parts[1], pairno, file_parts[3])
        if pairname in files:
          files.pop( files.index(pairname) )
          verified_files.append(file_parts[0])
          verified_files.append(pairname)
    if verified_files == []:
      self.logger.error("No correctly named fastq files found in directory. Exited.")
      sys.exit()
    return verified_files
 
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

  def create_trimjob(self):
    batchfile = open(self.batchfile, "a+")
    files = self.verify_fastq()
    i=0
    j=1
    while i < len(files):
      outfile = files[i].split('.')[0][:-2]
      if not outfile in self.trimmed_files:
        self.trimmed_files[outfile] = dict()
      self.trimmed_files[outfile]['fp'] = "{}/{}_trim_front_pair.fq".format(self.outdir, outfile)
      self.trimmed_files[outfile]['fu'] = "{}/{}_trim_front_unpair.fq".format(self.outdir, outfile)
      self.trimmed_files[outfile]['rp'] = "{}/{}_trim_rev_pair.fq".format(self.outdir, outfile)
      self.trimmed_files[outfile]['ru'] = "{}/{}_trim_rev_unpair.fq".format(self.outdir, outfile)
      
      batchfile.write("# Trimmomatic set {}\n".format(j))
      batchfile.write("trimmomatic-0.36.jar PE -threads {} -phred33 {}/{} {}/{} {} {} {} {}\
      ILLUMINACLIP:{}/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n\n"\
      .format(self.config["slurm_header"]["threads"], self.indir, files[i], self.indir, files[i+1],\
      self.trimmed_files[outfile]['fp'], self.trimmed_files[outfile]['fu'], self.trimmed_files[outfile]['rp'], self.trimmed_files[outfile]['ru'], self.config["folders"]["installations"]))
      i=i+2
      j+=1
    batchfile.close()

  def interlace_files(self):
    """Interlaces all unpaired files"""
    suffix = "_unpaired_interlaced.fq"
    for name, v in self.trimmed_files.items():
      self.logger.info("Creating unpaired interlace file for run {}".format(name))
      f = open("{}/{}{}".format(self.outdir, name, suffix), "w")
      pdb.set_trace()
      fu = open("{}".format(v['fu']), "r")
      ru = open("{}".format(v['ru']), "r")
      f.write(fu.read())
      f.write(ru.read())
      f.close()
      fu.close()
      ru.close()
      self.trimmed_files[name]['i'] = "{}/{}{}".format(self.outdir, name, suffix)
    
  def create_spadesjob(self):
    batchfile = open(self.batchfile, "a+")
    #memory is actually 128 per node regardless of cores.
    batchfile.write("# Spades assembly\n")
    batchfile.write("spades.py --threads {} --memory {} -o {}/assembly"\
    .format(self.config["slurm_header"]["threads"], 8*int(self.config["slurm_header"]["threads"]), self.outdir))
    
    libno = 1
    for k,v in self.trimmed_files.items():
      batchfile.write(" --pe{}-1 {}".format(libno, self.trimmed_files[k]['fp']))
      batchfile.write(" --pe{}-2 {}".format(libno, self.trimmed_files[k]['rp']))
      # Method requires manager
      #batchfile.write(" --pe{}-s {}".format(libno, self.trimmed_files[k]['i']))
      libno += 1

    batchfile.write("\n\n")
    batchfile.close()

  def create_blastjob(self):
    #index database
    batchfile = open(self.batchfile, "a+")

    #Establish organism
    refname = ""
    indexed = 0
    refs = os.listdir(self.config["folders"]["references"])
    for file in refs:
      hit = re.search('({}\w*).xmfa'.format(self.organism), file)
      if hit:
        refname = hit.group(1)
      if re.search('{}\w*.\w+'.format(self.organism), file):
        indexed = indexed + 1
    if refname == "":
      self.logger.error("Bad reference name given, no reference database found!")
      sys.exit()

    if indexed < 5:
      batchfile.write("# Blast database indexing. Only necessary for initial run of organism\n")
      batchfile.write("cd {} && makeblastdb -in {}/{}.xmfa -dbtype nucl -parse_seqids -out {}\n\n".format(\
      self.config["folders"]["references"], self.config["folders"]["references"], refname, refname))
    #create run
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send\""
    batchfile.write("# BLAST MLST alignment\n")
    batchfile.write("blastn -db {}/{} -query {}/assembly/contigs.fasta -out {}/loci_query_tab.txt -num_threads {} -max_target_seqs 1 -outfmt {}\n\n".format(\
    self.config["folders"]["references"], refname, self.outdir, self.outdir, self.config["slurm_header"]["threads"], blast_format))

#TODO: This function will be managed outside
  def create_job(self):
    if not os.path.exists(self.outdir):
      os.makedirs(self.outdir)
    self.batchfile = "{}/runfile.sbatch".format(self.outdir)

    self.create_header()
    self.create_trimjob()
    # Method requires manager
    #self.interlace_files()
    self.create_spadesjob()
    self.create_blastjob()
    self.logger.info("Created runfile for project {} in folder {}".format(self.name, self.outdir))
