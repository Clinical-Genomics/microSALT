"""Creates sbatch jobs for MLST instances
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import glob
import os
import re
import shutil
import subprocess
import time

from datetime import datetime
from microSALT.store.lims_fetcher import LIMS_Fetcher
from microSALT.store.db_manipulator import DB_Manipulator
#from microSALT.utils.scraper import Scraper
#from microSALT.utils.reporter import Reporter

class Job_Creator():

  def __init__(self, indir, config, log, finishdir="", timestamp=""):
    self.config = config
    self.logger = log
    self.batchfile = ""
    self.indir = os.path.abspath(indir)
    self.name = os.path.basename(os.path.normpath(indir))

    self.now = timestamp
    if timestamp != "":
      self.now = timestamp
      temp = timestamp.replace('_','.').split('.')
      self.dt = datetime(int(temp[0]),int(temp[1]),int(temp[2]),int(temp[3]),int(temp[4]),int(temp[5]))
    else:
      self.dt = datetime.now()
      self.now = time.strftime("{}.{}.{}_{}.{}.{}".\
      format(self.dt.year, self.dt.month, self.dt.day, self.dt.hour, self.dt.minute, self.dt.second))

    #Attempting writing on slurm
    self.outdir = "/scratch/$SLURM_JOB_ID/workdir/{}_{}".format(self.name, self.now)
    self.finishdir = finishdir
    if self.finishdir == "":
      self.finishdir="{}/{}_{}".format(config["folders"]["results"], self.name, self.now)

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
                 self.config["slurm_header"]["job_prefix"], self.name,self.config["slurm_header"]["qos"], self.finishdir, self.name)
    return headerline

  def verify_fastq(self):
    """ Uses arg indir to return a list of PE fastq tuples fulfilling naming convention """
    files = os.listdir(self.indir)
    if files == []:
      raise Exception("Directory {} lacks fastq files.".format(self.indir))
    verified_files = list()
    for file in files:
      file_match = re.match( self.config['regex']['file_pattern'], file)
      if file_match:
        #Check that symlinks resolve
        path = '{}/{}'.format(self.indir, file)
        if os.path.islink(path):
          if not os.path.exists(os.readlink(path)):
            raise Exception("Some fastq files are unresolved symlinks in directory {}.".format(self.indir))

        #Make sure both mates exist
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
          raise Exception("Some fastq files have no mate in directory {}.".format(self.indir))
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
    batchfile.write("spades.py --threads {} --careful --memory {} -o {}/assembly"\
    .format(self.config["slurm_header"]["threads"], 8*int(self.config["slurm_header"]["threads"]), self.outdir))

    libno = 1
    for k,v in self.trimmed_files.items():
      batchfile.write(" --pe{}-1 {}".format(libno, self.trimmed_files[k]['fp']))
      batchfile.write(" --pe{}-2 {}".format(libno, self.trimmed_files[k]['rp']))
      batchfile.write(" --pe{}-s {}".format(libno, self.trimmed_files[k]['i']))
      libno += 1

    batchfile.write("\n\n")
    batchfile.close()

  def create_resistancesection(self):
    """Creates a blast job for instances where many loci definition files make up an organism"""

    #Create run
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir {}/resistance\n\n".format(self.outdir))
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""
    res_list = glob.glob("{}/*.fsa".format(self.config["folders"]["resistances"]))
    for entry in res_list:
      batchfile.write("# BLAST Resistance search in {} for {}\n".format(self.organism, os.path.basename(entry[:-4])))
      batchfile.write("blastn -db {}  -query {}/assembly/contigs.fasta -out {}/resistance/{}.txt -task megablast -num_threads {} -outfmt {}\n".format(\
      entry[:-4], self.outdir, self.outdir, os.path.basename(entry[:-4]), self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()

  def create_blastsection(self):
    """Creates a blast job for instances where many loci definition files make up an organism"""

    #Create run
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir {}/blast\n\n".format(self.outdir))
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""
    tfa_list = glob.glob("{}/{}/*.tfa".format(self.config["folders"]["references"], self.organism))
    for entry in tfa_list:
      batchfile.write("# BLAST MLST alignment for {}, {}\n".format(self.organism, os.path.basename(entry[:-4])))
      batchfile.write("blastn -db {}  -query {}/assembly/contigs.fasta -out {}/blast/loci_query_{}.txt -task megablast -num_threads {} -outfmt {}\n".format(\
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
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir {}\n\n".format(trimdir))
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
    batchfile.write("mkdir {}/quast\n".format(self.outdir))
    batchfile.write("quast.py {}/assembly/contigs.fasta -o {}/quast\n\n".format(self.outdir, self.outdir))
    batchfile.close()

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
    proj_col['Customer_ID'] = self.lims_fetcher.data['Customer_ID']
    self.db_pusher.add_rec(proj_col, 'Projects')

  def create_sample(self, name):
    """Creates sample in database"""
    try:
      self.lims_fetcher.load_lims_sample_info(name)
      sample_col = self.db_pusher.get_columns('Samples')
      sample_col['CG_ID_sample'] = self.lims_fetcher.data['CG_ID_sample']
      sample_col['CG_ID_project'] = self.lims_fetcher.data['CG_ID_project']
      sample_col['Customer_ID_sample'] = self.lims_fetcher.data['Customer_ID_sample']
      sample_col["date_analysis"] = self.dt
      sample_col['organism']=self.lims_fetcher.data['organism']
      sample_col["application_tag"] = self.lims_fetcher.data['application_tag']
      sample_col["priority"] = self.lims_fetcher.data['priority']
      sample_col["date_sequencing"] = self.lims_fetcher.data['date_sequencing']
      sample_col["date_libprep"] = self.lims_fetcher.data['date_libprep']
      #self.db_pusher.purge_rec(sample_col['CG_ID_sample'], 'sample')
      self.db_pusher.add_rec(sample_col, 'Samples')
    except Exception as e:
      self.logger.error("Unable to add sample {} to database".format(self.name))

  def project_job(self, single_sample=False):
    if 'dry' in self.config and self.config['dry']==True:
      dry=True
    else:
      dry=False
    jobarray = list()
    if not os.path.exists(self.finishdir):
      os.makedirs(self.finishdir)
    try:
       if single_sample:
         self.create_project(os.path.normpath(self.indir).split('/')[-2])
       else:
        self.create_project(self.name)
    except Exception as e:
      self.logger.error("LIMS interaction failed. Unable to read/write project {}".format(self.name))
      #Start every sample job
    if single_sample:
      try:
        self.sample_job()
        headerargs = self.get_headerargs()
        outfile = self.get_sbatch()
        bash_cmd="sbatch {} {}".format(headerargs, outfile)
        if not dry and outfile != "":
          samproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
          output, error = samproc.communicate()
          jobno = re.search('(\d+)', str(output)).group(0)
          jobarray.append(jobno)
        else:
          self.logger.info("Suppressed command: {}".format(bash_cmd))
      except Exception as e:
        self.logger.error("Unable to analyze single sample {}".format(self.name))
    else:
      for (dirpath, dirnames, filenames) in os.walk(self.indir):
        for dir in dirnames:
          try:
            sample_in = "{}/{}".format(dirpath, dir)
            sample_out = "{}/{}".format(self.finishdir, dir)
            sample_instance = Job_Creator(sample_in, self.config, self.logger, sample_out, self.now)
            sample_instance.sample_job()
            headerargs = sample_instance.get_headerargs()
            outfile = ""
            if os.path.isfile(sample_instance.get_sbatch()):
              outfile = sample_instance.get_sbatch()
            bash_cmd="sbatch {} {}".format(headerargs, outfile)
            if not dry and outfile != "":
              projproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
              output, error = projproc.communicate()
              jobno = re.search('(\d+)', str(output)).group(0)
              jobarray.append(jobno)
            else:
              self.logger.info("Suppressed command: {}".format(bash_cmd))
          except Exception as e:
            pass
    if not dry:
      self.finish_job(jobarray, single_sample)

  def finish_job(self, joblist, single_sample=False):
    """ Uploads data and sends an email once all analysis jobs are complete.
        Sbatch > Srun to avoid issues with 50+ jobs """

    startfile = "{}/run_started.out".format(self.finishdir)
    mailfile = "{}/mailjob.sh".format(self.finishdir)
    mb = open(mailfile, "w+")
    sb = open(startfile, "w+")
    sb.write("#!/bin/sh\n\n")
    sb.close()
    mb.write("#!/bin/sh\n\n")
    mb.write("#Uploading of results to database and production of report\n")
    if 'MICROSALT_CONFIG' in os.environ:
      mb.write("export MICROSALT_CONFIG={}\n".format(os.environ['MICROSALT_CONFIG']))
    mb.write("source activate $CONDA_DEFAULT_ENV\n")
    if not single_sample:
      mb.write("microSALT finish project {} --input {} --rerun --email {}\n".\
               format(self.name, self.finishdir, self.config['regex']['mail_recipient']))
    else:
      mb.write("microSALT finish sample {} --input {} --rerun --email {}\n".\
               format(self.name, self.finishdir, self.config['regex']['mail_recipient']))
    mb.write("touch {}/run_complete.out".format(self.finishdir))
    mb.close()

    massagedJobs = list()
    final = ':'.join(joblist)
    #Create subtracker if more than 50 samples
    maxlen = 50
    if len(joblist) > maxlen:
      i = 1
      while i <= len(joblist):
        if i+maxlen < len(joblist):
          massagedJobs.append(':'.join(joblist[i-1:i+maxlen-1]))
        else:
          massagedJobs.append(':'.join(joblist[i-1:-1]))
        i += maxlen
      for entry in massagedJobs:
        if massagedJobs.index(entry) < len(massagedJobs)-1:
          head = "-A {} -p core -n 1 -t 00:00:10 -J {}_{}_SUBTRACKER --qos {} --dependency=afterany:{}"\
                 .format(self.config["slurm_header"]["project"],self.config["slurm_header"]["job_prefix"],\
                         self.name,self.config["slurm_header"]["qos"],entry)
          bash_cmd="sbatch {} {}".format(head, startfile)
          mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
          output, error = mailproc.communicate()
          jobno = re.search('(\d+)', str(output)).group(0)
          massagedJobs[massagedJobs.index(entry)+1] += ":{}".format(jobno)
        else:
          final = entry
          break

    head = "-A {} -p core -n 1 -t 06:00:00 -J {}_{}_MAILJOB --qos {} --open-mode append --dependency=afterany:{} --output {}"\
            .format(self.config["slurm_header"]["project"],self.config["slurm_header"]["job_prefix"],\
                    self.name,self.config["slurm_header"]["qos"],\
           final, self.config['folders']['log_file'],  self.config['regex']['mail_recipient'])
    bash_cmd="sbatch {} {}".format(head, mailfile)
    mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
    output, error = mailproc.communicate()

  def sample_job(self):
    try:
      self.trimmed_files = dict()
      if not os.path.exists(self.finishdir):
        os.makedirs(self.finishdir)
      try:
        self.organism = self.lims_fetcher.get_organism_refname(self.name, external=False)
        # This is one job
        self.batchfile = "{}/runfile.sbatch".format(self.finishdir)
        batchfile = open(self.batchfile, "w+")
        batchfile.write("#!/bin/sh\n\n")
        batchfile.write("mkdir -p {}\n".format(self.outdir))
        batchfile.close()

        self.create_trimsection()
        self.interlace_files()
        self.create_spadessection()
        self.create_quastsection()
        self.create_blastsection()
        self.create_resistancesection()
        batchfile = open(self.batchfile, "a+")
        batchfile.write("cp -r {}/* {}".format(self.outdir, self.finishdir))
        batchfile.close()

        self.logger.info("Created runfile for sample {} in folder {}".format(self.name, self.outdir))
      except Exception as e:
        raise
      try:
        self.create_sample(self.name)
      except Exception as e:
        self.logger.error("Unable to access LIMS info for sample {}".format(self.name))
    except Exception as e:
      self.logger.error("Unable to create job for sample {}\nSource: {}".format(self.name, str(e)))
      shutil.rmtree(self.finishdir, ignore_errors=True)
      raise
