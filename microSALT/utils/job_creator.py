"""Creates sbatch jobs for MLST instances
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import glob
import os
import re
import shutil
import subprocess
import sys
import time
import yaml

from microSALT.store.lims_fetcher import LIMS_Fetcher
from microSALT.store.db_manipulator import DB_Manipulator

class Job_Creator():

  def __init__(self, indir, config, log, outdir="", timestamp=""):
    self.config = config
    self.logger = log
    self.batchfile = ""
    self.indir = os.path.abspath(indir)
    self.name = os.path.basename(os.path.normpath(indir))

    self.now = timestamp
    if timestamp == "":
      self.now = time.strftime("%Y.%m.%d_%H.%M.%S")
    self.outdir = outdir
    if self.outdir == "":
      self.outdir="{}/{}_{}".format(config["folders"]["results"], os.path.basename(os.path.normpath(self.indir)), self.now)
    
    self.db_pusher=DB_Manipulator(config, log)
    self.trimmed_files = dict()
    self.organism = ""
    self.lims_fetcher = LIMS_Fetcher(config, log)

  def get_sbatch(self):
    """ Returns sbatchfile, slightly superflous"""
    return self.batchfile

  def get_headerargs(self):
    headerline = "-A {} -p {} -n {} -t {} -J {}_{} --qos {} --output {}/slurm_{}.log".format(self.config["slurm_header"]["project"],\
                 self.config["slurm_header"]["type"], self.config["slurm_header"]["threads"],self.config["slurm_header"]["time"],\
                 self.config["slurm_header"]["job_prefix"], self.name,self.config["slurm_header"]["qos"], self.outdir, self.name)
    return headerline

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
          pairname = "{}{}{}".format(file_match.string[:file_match.end(1)-1] , pairno, \
                      file_match.string[file_match.end(1):file_match.end()])
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
 
  def interlace_files(self):
    """Interlaces all unpaired files"""
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# Interlaced unpaired reads file creation\n")
    suffix = "_unpaired_interlaced.fq"
    for name, v in self.trimmed_files.items():
      interfile = "{}/trimmed/{}{}".format(self.outdir, name, suffix)
      # Spammed a bit too much
      #self.logger.info("Created unpaired interlace file for sample {}".format(name))
      batchfile.write("touch {}\n".format(interfile))
      batchfile.write("cat {} >> {}\n".format(v['fu'], interfile))
      batchfile.write("cat {} >> {}\n".format(v['ru'], interfile))
      self.trimmed_files[name]['i'] = "{}".format(interfile)
      batchfile.write("\n")
    batchfile.close()    

  def create_spadessection(self):
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

  def index_db(self, full_dir, suffix):
    """Check for indexation, makeblastdb job if not enough of them."""
    try:
      batchfile = open(self.batchfile, "a+")
      files = os.listdir(full_dir)
      sufx_files = glob.glob("{}/*{}".format(full_dir, suffix)) #List of source files
      nin_suff = sum([1 for elem in files if 'nin' in elem]) #Total number of nin files 
      if nin_suff < len(sufx_files):
        batchfile.write("# Blast database indexing. Only necessary for initial run against target\n")
        for file in sufx_files:
          if '.fsa' in suffix:
            batchfile.write("cd {} && makeblastdb -in {}/{} -dbtype nucl -out {}\n".format(\
            full_dir, full_dir, os.path.basename(file),  os.path.basename(file[:-4])))
          else:
            batchfile.write("cd {} && makeblastdb -in {}/{} -dbtype nucl -parse_seqids -out {}\n".format(\
            full_dir, full_dir, os.path.basename(file),  os.path.basename(file[:-4])))
          ref = open(file, "r")
          ref.readline()
      batchfile.write("\n")
      batchfile.close()
    except Exception as e:
      self.logger.error("Unable to index {} of organism {} against reference {}".format(self.name, self.organism, file))

  def create_cgmlstsection(self):
    """Creates a blast job against a known reference/expanded geneset"""
    self.index_db("{}".format(self.config["folders"]["gene_set"]), '.gst')
    if not os.path.exists("{}/cgmlst".format(self.outdir)):
      os.makedirs("{}/cgmlst".format(self.outdir))

    #Create run
    batchfile = open(self.batchfile, "a+")
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length qseq\""
    entry = "{}/{}.gst".format(self.config["folders"]["gene_set"], self.organism)
    
    batchfile.write("# BLAST cgMLST search in {} for {}\n".format(self.organism, os.path.basename(entry[:-4])))
    batchfile.write("blastn -db {} -query {}/assembly/contigs.fasta -out {}/cgmlst/loci_query_cgmlst.txt -task megablast -num_threads {} -perc_identity 95 -outfmt {}\n".format(\
    entry[:-4], self.outdir, self.outdir, self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()

  def create_resistancesection(self):
    """Creates a blast job for each resistance gene of an organism"""
    self.index_db("{}".format(self.config["folders"]["resistances"]), '.fsa')
    if not os.path.exists("{}/resistance".format(self.outdir)):
      os.makedirs("{}/resistance".format(self.outdir))

    #Create run
    batchfile = open(self.batchfile, "a+")
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""
    res_list = glob.glob("{}/*.fsa".format(self.config["folders"]["resistances"]))
    for entry in res_list:
      batchfile.write("# BLAST Resistance search in {} for {}\n".format(self.organism, os.path.basename(entry[:-4])))
      batchfile.write("blastn -db {}  -query {}/assembly/contigs.fasta -out {}/resistance/{}.txt -task megablast -num_threads {} -max_target_seqs 1 -outfmt {}\n".format(\
      entry[:-4], self.outdir, self.outdir, os.path.basename(entry[:-4]), self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()

  def create_blastsection(self):
    """Creates a blast job for instances where many loci definition files make up an organism"""
    self.index_db("{}/{}".format(self.config["folders"]["references"], self.organism), '.tfa')
    if not os.path.exists("{}/blast".format(self.outdir)):
      os.makedirs("{}/blast".format(self.outdir))
    
    #Create run
    batchfile = open(self.batchfile, "a+")
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""
    tfa_list = glob.glob("{}/{}/*.tfa".format(self.config["folders"]["references"], self.organism))
    for entry in tfa_list:
      batchfile.write("# BLAST MLST alignment for {}, {}\n".format(self.organism, os.path.basename(entry[:-4])))
      batchfile.write("blastn -db {}  -query {}/assembly/contigs.fasta -out {}/blast/loci_query_{}.txt -task megablast -num_threads {} -max_target_seqs 1 -outfmt {}\n".format(\
      entry[:-4], self.outdir, self.outdir, os.path.basename(entry[:-4]), self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()

  def create_trimsection(self):
    for root, dirs, files in os.walk(self.config["folders"]["adapters"]):
      if not "NexteraPE-PE.fa" in files:
        self.logger.error("Adapters folder at {} does not contain NexteraPE-PE.fa. Review paths.yml")
      else:
        break
    trimdir = "{}/trimmed".format(self.outdir)
    files = self.verify_fastq()
    if not os.path.exists(trimdir):
      os.makedirs(trimdir)
    batchfile = open(self.batchfile, "a+")
    i=0
    j=1
    while i < len(files):
      outfile = files[i].split('.')[0][:-2]
      if not outfile in self.trimmed_files:
        self.trimmed_files[outfile] = dict()
      self.trimmed_files[outfile]['fp'] = "{}/{}_trim_front_pair.fq".format(trimdir, outfile)
      self.trimmed_files[outfile]['fu'] = "{}/{}_trim_front_unpair.fq".format(trimdir, outfile)
      self.trimmed_files[outfile]['rp'] = "{}/{}_trim_rev_pair.fq".format(trimdir, outfile)
      self.trimmed_files[outfile]['ru'] = "{}/{}_trim_rev_unpair.fq".format(trimdir, outfile)

      batchfile.write("# Trimmomatic set {}\n".format(j))
      batchfile.write("trimmomatic PE -threads {} -phred33 {}/{} {}/{} {} {} {} {}\
      ILLUMINACLIP:{}NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\n\n"\
      .format(self.config["slurm_header"]["threads"], self.indir, files[i], self.indir, files[i+1],\
      self.trimmed_files[outfile]['fp'], self.trimmed_files[outfile]['fu'], self.trimmed_files[outfile]['rp'],\
      self.trimmed_files[outfile]['ru'], self.config["folders"]["adapters"]))
      i=i+2
      j+=1

  def create_quastsection(self):
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# QUAST QC metrics\n")
    batchfile.write("quast.py {}/assembly/contigs.fasta -o {}/quast\n".format(self.outdir, self.outdir))
    batchfile.write("\n")
    batchfile.close()
    if not os.path.exists("{}/quast".format(self.outdir)):
      os.makedirs("{}/quast".format(self.outdir))

  def create_project(self, name):
    """Creates project in database"""
    try:
      self.lims_fetcher.load_lims_project_info(name)
    except Exception as e:
      self.logger.error("Unable to load LIMS info for project {}".format(name))
    proj_col=dict()
    proj_col['CG_ID_project'] = name
    proj_col['Customer_ID_project'] = self.lims_fetcher.data['Customer_ID_project']
    proj_col['date_ordered'] = self.lims_fetcher.data['date_received']
    self.db_pusher.add_rec(proj_col, 'Projects')

  def create_sample(self, name):
    """Creates sample in database"""
    try:
      self.lims_fetcher.load_lims_sample_info(name)
      sample_col = self.db_pusher.get_columns('Samples') 
      sample_col['CG_ID_sample'] = self.lims_fetcher.data['CG_ID_sample']
      sample_col['CG_ID_project'] = self.lims_fetcher.data['CG_ID_project']
      sample_col['Customer_ID_sample'] = self.lims_fetcher.data['Customer_ID_sample']
      sample_col["date_analysis"] = self.now
      self.db_pusher.add_rec(sample_col, 'Samples')
    except Exception as e:
      self.logger.error("Unable to add sample {} to database".format(self.name))

  def project_job(self, single_sample=False):
    jobarray = list()
    if not os.path.exists(self.outdir):
      os.makedirs(self.outdir)

    try:
       if single_sample:
         self.create_project(os.path.normpath(self.indir).split('/')[-2])
       else:
        self.create_project(self.name)
    except Exception as e:
      self.logger.error("LIMS interaction failed. Unable to read/write project {}".format(self.name)) 
    try:
      #Start every sample job
      if single_sample:
        self.sample_job()
        headerargs = self.get_headerargs()
        outfile = self.get_sbatch()
        bash_cmd="sbatch {} {}".format(headerargs, outfile)
        process = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        jobno = re.search('(\d+)', str(output)).group(0)
        jobarray.append(jobno)
      else:
        for (dirpath, dirnames, filenames) in os.walk(self.indir):
          for dir in dirnames:
            sample_in = "{}/{}".format(dirpath, dir)
            sample_out = "{}/{}".format(self.outdir, dir)
            sample_instance = Job_Creator(sample_in, self.config, self.logger, sample_out, self.now) 
            sample_instance.sample_job()
            headerargs = sample_instance.get_headerargs()
            outfile = sample_instance.get_sbatch()
            bash_cmd="sbatch {} {}".format(headerargs, outfile)
            output, error = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE).communicate()
            jobno = re.search('(\d+)', str(output)).group(0)
            jobarray.append(jobno)
      #Mail job
      mailline = "srun -A {} -p core -n 1 -t 00:00:10 -J {}_{}_TRACKER --qos {} --dependency=afterany:{} --output {}/run_complete.out --mail-user={} --mail-type=END pwd"\
                 .format(self.config["slurm_header"]["project"],self.config["slurm_header"]["job_prefix"], self.name,self.config["slurm_header"]["qos"],\
                 ':'.join(jobarray), self.outdir,  self.config['regex']['mail_recipient'])
      mailarray = mailline.split()
      subprocess.Popen(mailarray, stdin=None, stdout=None, stderr=None)
    except Exception as e:
      self.logger.warning("Failed to spawn project at {}\nSource: {}".format(self.outdir, str(e)))
      shutil.rmtree(self.outdir, ignore_errors=True)

  def sample_job(self):
    self.trimmed_files = dict()
    if not os.path.exists(self.outdir):
      os.makedirs(self.outdir)
    try:
      self.organism = self.lims_fetcher.get_organism_refname(self.name, external=False)
      # This is one job 
      self.batchfile = "{}/runfile.sbatch".format(self.outdir)
      batchfile = open(self.batchfile, "w+")
      batchfile.write("#!/bin/sh\n\n")
      batchfile.close()

      self.create_trimsection()
      self.interlace_files()
      self.create_spadessection()
      self.create_quastsection()
      self.create_blastsection()
      self.create_resistancesection()
      self.create_cgmlstsection()
      self.logger.info("Created runfile for sample {} in folder {}".format(self.name, self.outdir))
    except Exception as e:
      raise Exception("Unable to create job for instance {}\nSource: {}".format(self.indir, str(e)))
    try: 
      self.create_sample(self.name)
    except Exception as e:
      self.logger.error("LIMS interaction failed. Unable to read/write sample {}".format(self.name))
