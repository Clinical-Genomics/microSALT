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
    self.concat_files = dict()
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
    """Interlaces all trimmed files"""
    fplist = list()
    kplist = list()
    ilist = list()
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# Interlaced trimmed files\n")
 
    for name, v in self.trimmed_files.items():
      fplist.append( v['fp'] )
      kplist.append( v['rp'] )
      ilist.append( v['fu'] )
      ilist.append( v['ru'] )

    if len(kplist) != len(fplist) or len(ilist)/2 != len(kplist):
      raise Exception("Uneven distribution of trimmed files. Invalid trimming step {}".format(name))
    self.concat_files['f'] = "{}/trimmed/{}{}".format(self.outdir,self.name, "_trim_front_pair.fq")
    self.concat_files['r'] = "{}/trimmed/{}{}".format(self.outdir,self.name, "_trim_rev_pair.fq")
    self.concat_files['i'] = "{}/trimmed/{}{}".format(self.outdir,self.name, "_trim_unpaired.fq")
    for k, v in self.concat_files.items():
      batchfile.write("touch {}\n".format(v))
    
    batchfile.write("cat {} >> {}\n".format(' '.join(fplist), self.concat_files['f']))
    batchfile.write("cat {} >> {}\n".format(' '.join(kplist), self.concat_files['r']))
    batchfile.write("cat {} >> {}\n".format(' '.join(ilist), self.concat_files['i']))
    batchfile.write("rm {} {} {}\n".format(' '.join(fplist), ' '.join(kplist), ' '.join(ilist)))
    batchfile.write("\n")
    batchfile.close()    

  def create_assemblysection(self):
    batchfile = open(self.batchfile, "a+")
    #memory is actually 128 per node regardless of cores.
    batchfile.write("# Spades assembly\n")
    batchfile.write("spades.py --threads {} --careful --memory {} -o {}/assembly"\
    .format(self.config["slurm_header"]["threads"], 8*int(self.config["slurm_header"]["threads"]), self.outdir))
    
    batchfile.write(" --pe1-1 {}".format(self.concat_files['f']))
    batchfile.write(" --pe1-2 {}".format(self.concat_files['r']))
    batchfile.write(" --pe1-s {}".format(self.concat_files['i']))

    batchfile.write("\n\n")
    batchfile.close()

  def create_resistancesection(self):
    """Creates a blast job for instances where many loci definition files make up an organism"""
    #self.index_db("{}".format(self.config["folders"]["resistances"]), '.fsa')

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

  def create_mlstsection(self):
    """Creates a blast job for instances where many loci definition files make up an organism"""
    
    #Create run
    batchfile = open(self.batchfile, "a+")
    batchfile.write("mkdir {}/blast\n\n".format(self.outdir))
    blast_format = "\"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length\""
    tfa_list = glob.glob("{}/{}/*.tfa".format(self.config["folders"]["references"], self.organism))
    for entry in tfa_list:
      batchfile.write("# BLAST MLST alignment for {}, {}\n".format(self.organism, os.path.basename(entry[:-4])))
      batchfile.write("blastn -db {}  -query {}/assembly/contigs.fasta -out {}/blast/loci_query_{}.txt -task megablast -num_threads {} -max_target_seqs 1 -outfmt {}\n".format(\
      entry[:-4], self.outdir, self.outdir, os.path.basename(entry[:-4]), self.config["slurm_header"]["threads"], blast_format))
    batchfile.write("\n")
    batchfile.close()

  def create_variantsection(self):
    """ Creates a job for variant calling based on local alignment """
    ref = "{}/{}.fasta".format(self.config['folders']['genomes'],self.lims_fetcher.data['reference'])
    outbase = "{}/alignment/{}_{}".format(self.outdir, self.name, self.lims_fetcher.data['reference'])

    #Create run
    batchfile = open(self.batchfile, "a+")
    batchfile.write("# Variant calling based on local alignment\n")
    batchfile.write("mkdir {}/alignment\n".format(self.outdir))
    batchfile.write("bwa mem -t {} -o {}.sam {} {} {}\n".format(self.config["slurm_header"]["threads"], outbase, ref ,self.concat_files['f'], self.concat_files['r']))
    batchfile.write("samtools view --threads {} -b -o {}.bam -T {} {}.sam\n".format(self.config["slurm_header"]["threads"], outbase, ref, outbase))
    batchfile.write("samtools sort --threads {} -n -o {}.bam_sort {}.bam\n".format(self.config["slurm_header"]["threads"], outbase, outbase))
    batchfile.write("samtools fixmate --threads {} -m {}.bam_sort {}.bam_sort_ms\n".format(self.config["slurm_header"]["threads"], outbase, outbase))
    batchfile.write("samtools sort --threads {} -n -o {}.bam_sort {}.bam_sort_ms\n".format(self.config["slurm_header"]["threads"], outbase, outbase))
    batchfile.write("samtools markdup -r -s --threads {} --reference {} --output-fmt bam {}.bam_sort {}.bam_sort_mkdup\n".format(self.config["slurm_header"]["threads"], ref, outbase, outbase))
    batchfile.write("samtools rmdup --reference {} {}.bam_sort_mkdup {}.bam_sort_rmdup\n".format(ref, outbase, outbase))
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

  def create_assemblystats_section(self):
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
    try:
      #Start every sample job
      if single_sample:
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
      else:
        for (dirpath, dirnames, filenames) in os.walk(self.indir):
          for dir in dirnames:
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
      if not dry:
        self.finish_job(jobarray)
    except Exception as e:
      self.logger.error("Issues handling some samples of project at {}\nSource: {}".format(self.finishdir, str(e)))
      #shutil.rmtree(self.finishdir, ignore_errors=True)

  def finish_job(self, joblist):
    """ Uploads data and sends an email once all analysis jobs are complete. """

    startfile = "{}/run_started.out".format(self.finishdir)
    mailfile = "{}/mailjob.sh".format(self.finishdir)
    mb = open(mailfile, "w+")
    sb = open(startfile, "w+")
    sb.write("#!/bin/sh\n\n")
    sb.close()
    mb.write("#!/bin/sh\n\n")
    mb.write("#Uploading of results to database and production of report\n")
    #mb.write("source ~/.bash_profile\n")
    if 'MICROSALT_CONFIG' in os.environ:
      mb.write("export MICROSALT_CONFIG={}\n".format(os.environ['MICROSALT_CONFIG']))
    mb.write("source activate $CONDA_DEFAULT_ENV\n")
    if not len(joblist) == 1:
      mb.write("microSALT finish project {} --input {} --rerun\n".format(self.name, self.finishdir))
    else:
      mb.write("microSALT finish sample {} --input {} --rerun\n".format(self.name, self.finishdir))
    mb.write("Analysis done!\n")
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
                 .format(self.config["slurm_header"]["project"],self.config["slurm_header"]["job_prefix"], self.name,self.config["slurm_header"]["qos"],\
                 entry)
          bash_cmd="sbatch {} {}".format(head, startfile)
          mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
          output, error = mailproc.communicate()
          jobno = re.search('(\d+)', str(output)).group(0)
          massagedJobs[massagedJobs.index(entry)+1] += ":{}".format(jobno)
        else:
          final = entry
          break 

    head = "-A {} -p core -n 1 -t 06:00:00 -J {}_{}_MAILJOB --qos {} --dependency=afterany:{} --output {}/run_complete.out"\
            .format(self.config["slurm_header"]["project"],self.config["slurm_header"]["job_prefix"], self.name,self.config["slurm_header"]["qos"],\
           final, self.finishdir,  self.config['regex']['mail_recipient'])
    bash_cmd="sbatch {} {}".format(head, mailfile)
    mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
    output, error = mailproc.communicate()

  def sample_job(self):
    """ Writes necessary sbatch job for each individual sample """
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
      self.create_variantsection()
      self.create_assemblysection()
      self.create_assemblystats_section()
      self.create_mlstsection()
      self.create_resistancesection()
      batchfile = open(self.batchfile, "a+")
      batchfile.write("cp -r {}/* {}".format(self.outdir, self.finishdir))
      batchfile.close()

      self.logger.info("Created runfile for sample {} in folder {}".format(self.name, self.outdir))
    except Exception as e:
      self.logger.error("Unable to create job for instance {}\nSource: {}".format(self.indir, str(e)))
      shutil.rmtree(self.finishdir, ignore_errors=True) 
    try: 
      self.create_sample(self.name)
    except Exception as e:
      self.logger.error("Unable to access LIMS info for sample {}".format(self.name))
