"""Scrapes output files for data and adds them to the database
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import glob
import os
import re
import sys
import time
import yaml

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.store.lims_fetcher import LIMS_Fetcher
from microSALT.utils.reporter import Reporter
from microSALT.utils.job_creator import Job_Creator

# TODO: Rewrite so samples use seperate objects
class Scraper():

  def __init__(self, infolder, config, log):
    self.config = config
    self.logger = log
    self.infolder = os.path.abspath(infolder)
    self.sampledir = ""
   
    last_folder = os.path.basename(os.path.normpath(self.infolder)) 
    self.name = last_folder.split('_')[0]
    #TODO: Replace date from dir with entry from analysis files/database
    self.date = "{} {}".format(re.sub('\.','-', last_folder.split('_')[1]), re.sub('\.',':', last_folder.split('_')[2]))
    self.db_pusher=DB_Manipulator(config, log)
    self.lims_fetcher=LIMS_Fetcher(config, log)
    self.job_fallback=Job_Creator("", config, log)
    self.lims_sample_info = {}

  def scrape_project(self):
    """Scrapes a project folder for information"""
    if self.config['rerun']:
      self.db_pusher.purge_rec(self.name, 'project')
    if not self.db_pusher.exists('Projects', {'CG_ID_project':self.name}):
      self.logger.error("Re-filling project {}".format(self.name))
      self.job_fallback.create_project(self.name)

    #Scrape order matters a lot!
    for dir in os.listdir(self.infolder):
     if os.path.isdir("{}/{}".format(self.infolder, dir)): 
       self.sampledir = "{}/{}".format(self.infolder, dir)
       self.name = dir
       if not self.db_pusher.exists('Samples', {'CG_ID_sample':self.name}):
         self.logger.error("Re-filling sample {}".format(self.name))
         self.job_fallback.create_sample(self.name)
       self.scrape_all_loci()
       self.scrape_resistances()
       self.scrape_quast()

  def scrape_sample(self):
    """Scrapes a sample folder for information"""
    if self.config['rerun']:
      self.db_pusher.purge_rec(self.name, 'sample')

    self.lims_fetcher.load_lims_sample_info(self.name)
    if not self.db_pusher.exists('Projects', {'CG_ID_project':self.lims_fetcher.data['CG_ID_project']}):
      self.logger.error("Re-filling project {}".format(self.lims_fetcher.data['CG_ID_project']))
      self.job_fallback.create_project(self.lims_fetcher.data['CG_ID_project'])

    if not self.db_pusher.exists('Samples', {'CG_ID_sample':self.name}):
      self.logger.error("Re-filling sample {}".format(self.name))
      self.job_fallback.create_sample(self.name)

    #Scrape order matters a lot!
    self.sampledir = self.infolder
    self.scrape_all_loci()
    self.scrape_resistances()
    self.scrape_quast()

  def scrape_quast(self):
    """Scrapes a quast report for assembly information"""
    quast = dict()
    report = "{}/quast/report.tsv".format(self.sampledir)
    try:
      with open(report, 'r') as infile:
        for line in infile:
          lsplit = line.rstrip().split('\t')
          if lsplit[0] == '# contigs':
            quast['contigs'] = int(lsplit[1])
          elif lsplit[0] == 'Total length':
            quast['genome_length'] = int(lsplit[1])
          elif lsplit[0] == 'GC (%)':
            quast['gc_percentage'] = float(lsplit[1])
          elif lsplit[0] == 'N50':
            quast['n50'] = int(lsplit[1])

      self.db_pusher.upd_rec({'CG_ID_sample' : self.name}, 'Samples', quast)
      # Too spammy
      #self.logger.info("Project {} recieved quast stats: {}"\
      #               .format(self.name, quast))
    except Exception as e:
      self.logger.warning("Cannot generate quast statistics for {}".format(self.name))

  def get_locilength(self, analysis, reference, target):
    alleles=dict()
    targetPre = ">{}".format(target)
    lastallele=""
    if analysis=="Resistances":
      filename="{}/{}.fsa".format(self.config["folders"]["resistances"], reference)
    elif analysis=="Seq_types":
      loci = target.split('_')[0]
      filename="{}/{}/{}.tfa".format(self.config["folders"]["references"], reference, loci)

    f = open(filename,"r")
    for row in f:
      if ">" in row:
        lastallele = row.strip()
        alleles[lastallele] = ""
      else:
        alleles[lastallele] = alleles[lastallele] + row.strip()
    f.close()
    return len(alleles[targetPre])

  def scrape_resistances(self):
    q_list = glob.glob("{}/resistance/*".format(self.sampledir))
    for file in q_list:
      res_col = self.db_pusher.get_columns('Resistances')
      with open("{}".format(file), 'r') as sample:
        res_col["CG_ID_sample"] = self.name
        for line in sample:
          #Ignore commented fields
          if not line[0] == '#':
            elem_list = line.rstrip().split("\t")
            if not elem_list[1] == 'N/A':
              res_col["identity"] = elem_list[4]
              res_col["evalue"] = elem_list[5]
              res_col["bitscore"] = elem_list[6]
              res_col["contig_start"] = elem_list[7]
              res_col["contig_end"] = elem_list[8]
              res_col["subject_length"] =  elem_list[11]

              # Split elem 3 in loci (name) and allele (number) 
              res_col["gene"] = elem_list[3]
              res_col["instance"] = os.path.basename(file[:-4])
              res_col["span"] = float(res_col["subject_length"])/self.get_locilength('Resistances', res_col["instance"], res_col["gene"])

              # split elem 2 into contig node_NO, length, cov
              nodeinfo = elem_list[2].split('_')
              res_col["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
              res_col["contig_length"] = nodeinfo[3]
              res_col["contig_coverage"] = nodeinfo[5]
              self.db_pusher.add_rec(res_col, 'Resistances')

  def form_fingerprint(self, type):
    self.logger.info("Start cgMLST fingerprint generation. ETA 1-4 minutes per sample")
    if type == 'sample':
      self.lims_fetcher.load_lims_sample_info(self.name)
      self.sampledir = self.infolder
      self.scrape_cgmlst()
    elif type == 'project':
      for dir in os.listdir(self.infolder):
       if os.path.isdir("{}/{}".format(self.infolder, dir)):
         self.sampledir = "{}/{}".format(self.infolder, dir)
         self.name = dir
         start = time.time()
         self.scrape_cgmlst()
         self.logger.info("Time used for fingerprinting {}: {} seconds".format(self.name, int(time.time()-start)))
    else:
     self.logger.error("Invalid type programmed for cgMLST analysis")
     sys.exit(-1)

  def init_cgmlst(self):
    """Add downloaded gene list as default to database (if not already existing)"""
    organism = self.lims_fetcher.get_organism_refname(self.name)
    gene_list = "{}/{}.fna".format(self.config['folders']['gene_set'], organism)
    with open(gene_list, 'r') as gfile:
      for title in gfile:
        if title[0] == ">":
          title = title.replace('[','')
          title = title.replace(']','')
          title = title.replace('=',' ')
          #Title seperated
          title = title.split(' ')
          id = title[title.index('protein_id')+1]
          existing = self.db_pusher.query_rec('Profile_cgmlst', {'protein_id':id, 'organism':organism})
          if len(existing) < 1:
            self.db_pusher.add_rec({'protein_id':id, 'organism':organism, 'sequence':"Unassigned", 'allele':0}, 'Profile_cgmlst')

  def scrape_cgmlst(self):
    """Uploads cgmlst results to database and generates a fingerprint"""

    # self.init_cgmlst() <-- Unnecessary, only used for debugging. Even so needs revision, good idea though
    #Create sample specific node lookup table
    node_lookup = dict()
    ref_file = "{}/assembly/contigs.fasta".format(self.sampledir)
    with open (ref_file, 'r') as infile:
      for line in infile:
        if line[0] == ">":
          key = line[1:].rstrip()
          value = ""
        else:
          value += line.rstrip()
          node_lookup[key] = value

    #Associate each hit with it's sequence. BLAST displayed has some issues.
    allele_sequence = dict()
    res_file = "{}/cgmlst/loci_query_cgmlst.txt".format(self.sampledir)
    with open(res_file, 'r') as infile:
      for line in infile:
        if line[0] == "[":
          #Blast output seperated
          la = line.split('\t')
          title = la[0]
          title = title.replace('[','')
          title = title.replace(']','')
          title = title.replace('=',' ')
          #Title seperated
          title = title.split(' ')
          id = title[title.index('protein_id')+1]
          #Adds id entry with highest e-value
          if not id in allele_sequence:
            #la[-1] contains gaps
            if int(la[7]) < int(la[8]):
              #Indexes might be 1 away. Be careful
              allele_sequence[id] = node_lookup[la[2]][int(la[7]):int(la[8])]
            else:
              self.logger.error("Unhandled action. Index for {} come in reverse order".format(id))


    #Set all existing loci entries to N/A
    fp = dict()
    organism = self.lims_fetcher.get_organism_refname(self.name)

    #Check if found result exists in database
    for k,v in allele_sequence.items():
      push_dict = self.db_pusher.get_columns('Profile_cgmlst')
      top_index = self.db_pusher.top_index('Profile_cgmlst', {'protein_id':k, 'organism':organism}, 'allele')

      if top_index == -1:
        self.db_pusher.add_rec({'protein_id':k, 'organism':organism, 'sequence':v, 'allele':1}, 'Profile_cgmlst')
        fp[k] = 1
      else:
        exists = False
        for entry in self.db_pusher.query_rec('Profile_cgmlst', {'protein_id':k, 'organism':organism}):
          if entry.sequence in v or v in entry.sequence:
            fp[k] = entry.allele
            exists = True
            break
        if not exists:
          self.db_pusher.add_rec({'protein_id':k, 'organism':organism, 'sequence':v, 'allele':top_index+1}, 'Profile_cgmlst')
          fp[k] = top_index+1

    #Write fingerprint file
    fingerprint = "{}/cgmlst/fingerprint.txt".format(self.sampledir)
    fp_pointer = open(fingerprint,'w')
    for k, v in fp.items():
      fp_pointer.write("{}:{}\n".format(k,v))
    self.logger.info("Finished generating cgmlst results for {}".format(self.name)) 
  
  def scrape_all_loci(self):
    """Scrapes all BLAST output in a folder"""
    q_list = glob.glob("{}/blast/loci_query_*".format(self.sampledir))
    organism = self.lims_fetcher.get_organism_refname(self.name)
    self.db_pusher.upd_rec({'CG_ID_sample' : self.name}, 'Samples', {'organism': organism})

    for file in q_list:
      self.scrape_single_loci(file)
    #Requires all loci results to be initialized
    try:
      ST = self.db_pusher.alleles2st(self.name)
      self.db_pusher.upd_rec({'CG_ID_sample':self.name}, 'Samples', {'ST':ST})
      self.logger.info("Sample {} received ST {}".format(self.name, ST))
    except Exception as e:
      self.logger.error("{}".format(str(e)))

  def scrape_single_loci(self, infile):
    """Scrapes a single blast output file for MLST results"""
    seq_col = self.db_pusher.get_columns('Seq_types') 
    organism = self.lims_fetcher.get_organism_refname(self.name)
    if not os.path.exists(self.sampledir):
      self.logger.error("Invalid file path to infolder, {}".format(self.sampledir))
      sys.exit()
    try:
      with open("{}".format(infile), 'r') as insample:
        seq_col["CG_ID_sample"] = self.name

        for line in insample:
          #Ignore commented fields
          if not line[0] == '#':
            elem_list = line.rstrip().split("\t")
            if not elem_list[1] == 'N/A':
              seq_col["identity"] = elem_list[4]
              seq_col["evalue"] = elem_list[5]
              seq_col["bitscore"] = elem_list[6]
              seq_col["contig_start"] = elem_list[7]
              seq_col["contig_end"] = elem_list[8]
              seq_col["subject_length"] =  elem_list[11]
         
              # Split elem 3 in loci (name) and allele (number) 
              seq_col["loci"] = elem_list[3].split('_')[0]
              seq_col["allele"] = int(elem_list[3].split('_')[1])
              seq_col["span"] = float(seq_col["subject_length"])/self.get_locilength('Seq_types', organism, "{}_{}".format(seq_col["loci"], seq_col["allele"]))

              # split elem 2 into contig node_NO, length, cov
              nodeinfo = elem_list[2].split('_')
              seq_col["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
              seq_col["contig_length"] = nodeinfo[3]
              seq_col["contig_coverage"] = nodeinfo[5]
              self.db_pusher.add_rec(seq_col, 'Seq_types')
      # Too spammy
      #self.logger.info("Added allele {}={} of sample {} to table Seq_types".format(seq_col["loci"], seq_col["allele"], self.name))
    except Exception as e:
      self.logger.error("{}".format(str(e)))
