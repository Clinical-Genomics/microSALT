"""Creates sbatch jobs for MLST instances
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import glob
import os
import re
import shutil
import sys
import time
import yaml

from microSALT.store.lims_fetcher import LIMS_Fetcher

class Job_Creator():

  def __init__(self, indir, config, log, outdir=""):
    self.config = config
    self.logger = log
    self.now = time.strftime("%Y.%m.%d_%H.%M.%S")
    self.lims_fetcher = lims_fetcher=LIMS_Fetcher(config, log)
    if outdir == "":
      self.outdir="{}/{}_{}".format(config["folders"]["results"], os.path.basename(os.path.normpath(indir)), self.now)
    self.indir = os.path.abspath(indir)

    self.trimmed_files = dict()
    self.batchfile = ""
    self.organism = ""
    self.sample_name = os.path.basename(os.path.normpath(indir))

  def verify_fastq(self):
    """ Uses arg indir to return a list of PE fastq tuples fulfilling naming convention """
    files = os.listdir(self.indir)
    if files == []:
      raise Exception("Directory {} lacks fastq files.".format(self.indir))
    verified_files = list()
    for file in files:
      file_match = re.match( self.config['regex']['file_pattern'], file)
      #file_match = self.fileformat.match( file )
      if file_match:
        if file_match[1] == '1':
          pairno = '2'
          #Construct mate name
          pairname = "{}{}{}".format(file_match.string[:file_match.end(1)-1] , pairno, file_match.string[file_match.end(1):file_match.end()])
          if pairname in files:
            files.pop( files.index(pairname) )
            verified_files.append(file_match[0])
            verified_files.append(pairname)
        elif file_match[1] == '2':
          pass
        else:
          raise Exception("Some fastq files in directory have no mate in directory {}.".format(self.indir))
    if verified_files == []:
      raise Exception("No files in directory {} match file_pattern '{}'.".format(self.indir, self.config['regex']['file_pattern']))
    return verified_files
 
  def create_header(self):
    batchfile = open(self.batchfile, "w+")
    batchfile.write("#!/bin/bash -l\n\n")
    batchfile.write("#SBATCH -A {}\n".format(self.config["slurm_header"]["project"]))
    batchfile.write("#SBATCH -p {}\n".format(self.config["slurm_header"]["type"]))
    batchfile.write("#SBATCH -n {}\n".format(self.config["slurm_header"]["threads"]))
    batchfile.write("#SBATCH -t {}\n".format(self.config["slurm_header"]["time"]))
    batchfile.write("#SBATCH -J {}_{}\n".format(self.config["slurm_header"]["job_prefix"], self.sample_name))
    batchfile.write("#SBATCH --qos {}\n\n".format(self.config["slurm_header"]["qos"]))
    batchfile.close()

  def create_trimjob(self):
    for root, dirs, files in os.walk(self.config["folders"]["adapters"]):
      if not "NexteraPE-PE.fa" in files: 
        self.logger.error("Adapters folder at {} does not contain NexteraPE-PE.fa. Review paths.yml")
      else:
        break

    files = self.verify_fastq()
    batchfile = open(self.batchfile, "a+")
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
      batchfile.write("trimmomatic PE -threads {} -phred33 {}/{} {}/{} {} {} {} {}\
      ILLUMINACLIP:{}NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n\n"\
      .format(self.config["slurm_header"]["threads"], self.indir, files[i], self.indir, files[i+1],\
      self.trimmed_files[outfile]['fp'], self.trimmed_files[outfile]['fu'], self.trimmed_files[outfile]['rp'], self.trimmed_files[outfile]['ru'], self.config["folders"]["adapters"]))
      i=i+2
      j+=1
    batchfile.close()

  def interlace_files(self):
    """Interlaces all unpaired files"""
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# Interlaced unpaired reads file creation\n")
    suffix = "_unpaired_interlaced.fq"
    for name, v in self.trimmed_files.items():
      interfile = "{}/{}{}".format(self.outdir, name, suffix)
      self.logger.info("Creating unpaired interlace file for run {}".format(name))
      batchfile.write("touch {}\n".format(interfile))
      batchfile.write("cat {} >> {}\n".format(v['fu'], interfile))
      batchfile.write("cat {} >> {}\n".format(v['ru'], interfile))
      self.trimmed_files[name]['i'] = "{}".format(interfile)
      batchfile.write("\n")
    batchfile.close()    

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
      batchfile.write(" --pe{}-s {}".format(libno, self.trimmed_files[k]['i']))
      libno += 1

    batchfile.write("\n\n")
    batchfile.close()

  def index_db(self, full_dir):
    """Check for indexation, makeblastdb job if not enough of them."""
    try:
      batchfile = open(self.batchfile, "a+")
      files = os.listdir(full_dir)
      tfa_list = glob.glob("{}/*.tfa".format(full_dir))
      nin_suff = sum([1 for elem in files if 'nin' in elem]) #one type of index file 
      if nin_suff < len(tfa_list):
        batchfile.write("# Blast database indexing. Only necessary for initial run of organism\n")
        for file in tfa_list:
          batchfile.write("cd {} && makeblastdb -in {}/{} -dbtype nucl -parse_seqids -out {}\n".format(\
          full_dir, full_dir, os.path.basename(file),  os.path.basename(file[:-4])))
      batchfile.write("\n")
      batchfile.close()
    except Exception as e:
      self.logger.error("No associated reference for {} for specified organism: {}".format(self.sample_name, self.organism))

  def create_blastjob_single(self):
    """ Creates a blast job for instances where the definitions file is one per organism"""

    self.index_db(self.config["folders"]["references"])

    #create run
    batchfile = open(self.batchfile, "a+")
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send\""
    batchfile.write("# BLAST MLST alignment\n")
    batchfile.write("blastn -db {}/{} -query {}/assembly/contigs.fasta -out {}/loci_query_tab.txt -task megablast -num_threads {} -max_target_seqs 1 -outfmt {}\n\n".format(\
    self.config["folders"]["references"], self.organism, self.outdir, self.outdir, self.config["slurm_header"]["threads"], blast_format))
    batchfile.close()

  def create_blastjob_multi(self):
    """Creates a blast job for instances where many loci definition files make up an organism"""
    self.index_db("{}/{}".format(self.config["folders"]["references"], self.organism))
 
    #Create run
    batchfile = open(self.batchfile, "a+")
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send\""
    tfa_list = glob.glob("{}/{}/*.tfa".format(self.config["folders"]["references"], self.organism))
    for entry in tfa_list:
      batchfile.write("# BLAST MLST alignment for {}, {}\n".format(self.organism, os.path.basename(entry[:-4])))
      batchfile.write("blastn -db {}  -query {}/assembly/contigs.fasta -out {}/loci_query_{}.txt -task megablast -num_threads {} -max_target_seqs 1 -outfmt {}\n".format(\
      entry[:-4], self.outdir, self.outdir, os.path.basename(entry[:-4]), self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()


  def get_sbatch(self):
    return self.batchfile

  #TODO: Let project job spawn more job_creator objects rather than reassigning instance variables
  def project_job(self):
    proj_path = self.outdir
    if not os.path.exists(self.outdir):
      os.makedirs(self.outdir)
    concat_file = "{}/concatinated.sbatch".format(self.outdir)
    concat = open(concat_file, 'w+')
    try:
      concat.write("#!/bin/sh\n\n")
      for (dirpath, dirnames, filenames) in os.walk(self.indir):
        for dir in dirnames:
          self.outdir = "{}/{}".format(proj_path, dir) 
          self.indir = "{}/{}".format(dirpath, dir) 
          self.sample_name = dir 
          self.sample_job()
          outfile = self.get_sbatch()
          concat.write("sbatch {}\n".format(outfile))
      concat.close()
    except Exception as e:
      concat.close()
      self.logger.warning("Failed to spawn project at {}\nSource: {}".format(proj_path, str(e)))
      shutil.rmtree(proj_path, ignore_errors=True)

  def sample_job(self):
    self.trimmed_files = dict()
    try:
      self.organism = self.lims_fetcher.get_organism_refname(self.sample_name)
      if not os.path.exists(self.outdir):
        os.makedirs(self.outdir)
      self.batchfile = "{}/runfile.sbatch".format(self.outdir)
      
      self.create_header()
      self.create_trimjob()
      self.interlace_files()
      self.create_spadesjob()
      self.create_blastjob_multi()
      self.logger.info("Created runfile for project {} in folder {}".format(self.indir, self.outdir))
    except Exception as e:
      raise Exception("Unable to create job for instance {}\nSource: {}".format(self.indir, str(e)))
      #shutil.rmtree(self.outdir, ignore_errors=True)
