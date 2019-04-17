"""Scrapes output files for data and adds them to the database
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import glob
import os
import re
import string
import sys
import time

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
    self.gene2resistance = self.load_resistances()

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
       self.scrape_alignment()
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
    self.scrape_alignment()
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
    """ Find length of target within reference """
    alleles=dict()
    targetPre = ">{}".format(target)
    lastallele=""
    if analysis=="Resistances":
      filename="{}/{}.fsa".format(self.config["folders"]["resistances"], reference)
    elif analysis=="Seq_types":
      target = re.search('(.+)_(\w+)', target).group(1)
      filename="{}/{}/{}.tfa".format(self.config["folders"]["references"], reference, target)

    #Create dict with full name as key, associated nucleotides as value. 
    f = open(filename,"r")
    for row in f:
      if ">" in row:
        lastallele = row.strip()
        alleles[lastallele] = ""
      else:
        alleles[lastallele] = alleles[lastallele] + row.strip()
    f.close()
    try:
      return len(alleles[targetPre])
    except KeyError as e:
      self.logger.error("Target '{}' has been removed from current version of resFinder! Defaulting hit to length 1".format(targetPre))
      return 1

  def scrape_resistances(self):
    q_list = glob.glob("{}/resistance/*".format(self.sampledir))
    hypo = list()
    res_cols = self.db_pusher.get_columns('Resistances')
    for file in q_list:
      with open("{}".format(file), 'r') as sample:
        for line in sample:
          #Ignore commented fields
          if not line[0] == '#':
            elem_list = line.rstrip().split("\t")
            if not elem_list[1] == 'N/A':
              hypo.append(dict())
              hypo[-1]["CG_ID_sample"] = self.name
              hypo[-1]["identity"] = elem_list[4]
              hypo[-1]["evalue"] = elem_list[5]
              hypo[-1]["bitscore"] = elem_list[6]
              if int(elem_list[7]) < int(elem_list[8]):
                hypo[-1]["contig_start"] = int(elem_list[7])
                hypo[-1]["contig_end"] = int(elem_list[8])
              else:
                hypo[-1]["contig_start"] = int(elem_list[8])
                hypo[-1]["contig_end"] = int(elem_list[7])
              hypo[-1]["subject_length"] =  int(elem_list[11])
              hypo[-1]["instance"] = os.path.basename(file[:-4])

              # Split elem 3 in loci (name) and allele (number)
              partials = re.search('(.+)_(\d+){1,3}(?:_(\w+))*', elem_list[3])
              #if partials is None:
              #  print("DEBUG!")
              #  import pdb; pdb.set_trace()
              hypo[-1]["reference"] = partials.group(3)

              hypo[-1]["gene"] = partials.group(1)
              if hypo[-1]["gene"] in self.gene2resistance.keys():
                hypo[-1]["resistance"] = self.gene2resistance[hypo[-1]["gene"]]
              else:
                hypo[-1]["resistance"] = hypo[-1]["instance"].capitalize()

              hypo[-1]["span"] = float(hypo[-1]["subject_length"])/self.get_locilength('Resistances', hypo[-1]["instance"], partials.group(0))

              # split elem 2 into contig node_NO, length, cov
              nodeinfo = elem_list[2].split('_')
              hypo[-1]["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
              hypo[-1]["contig_length"] = int(nodeinfo[3])
              hypo[-1]["contig_coverage"] = nodeinfo[5]

    self.logger.info("{} candidate resistance hits found".format(len(hypo)))

    #Cleanup of overlapping hits
    ind = 0
    while ind < len(hypo)-1:
      targ = ind+1
      while targ < len(hypo):
        ignore = False
        if hypo[ind]["contig_name"] == hypo[targ]["contig_name"]:      
          #Overlapping or shared gene 
          if (hypo[ind]["contig_start"] >= hypo[targ]["contig_start"] and hypo[ind]["contig_start"] <= hypo[targ]["contig_end"]) or\
              (hypo[ind]["contig_end"] >= hypo[targ]["contig_start"] and hypo[ind]["contig_end"] <= hypo[targ]["contig_end"]) or \
              (hypo[ind]['gene'] == hypo[targ]['gene']):
            #Rightmost is worse
            if float(hypo[ind]["identity"])*(1-abs(1-hypo[ind]["span"])) >= float(hypo[targ]["identity"])*(1-abs(1-hypo[targ]["span"])):
              del hypo[targ]
              ignore = True
            #Leftmost is worse
            elif float(hypo[ind]["identity"])*(1-abs(1-hypo[ind]['span'])) < float(hypo[targ]["identity"])*(1-abs(1-hypo[targ]['span'])):
              del hypo[ind]
              targ = ind +1
              ignore = True
        if not ignore:
          targ += 1
        else:
          pass
      ind += 1

    self.logger.info("{} resistance hits were added after removing overlaps and duplicate hits".format( len(hypo)))
    for hit in hypo:
      self.db_pusher.add_rec(hit, 'Resistances')

  def load_resistances(self):
    """Loads common resistance names for genes"""
    conversions = dict()
    with open("{}/notes.txt".format(self.config["folders"]["resistances"])) as fh:
      for line in fh:
        if '#' not in line:
          line = line.split(':')
          cropped = re.sub(' resistance', '', line[1])
          conversions[line[0]] = cropped
          #Workaround for case issues
          conversions[line[0].lower()] = cropped
    return conversions

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
    hypo = list()
    seq_col = self.db_pusher.get_columns('Seq_types') 
    organism = self.lims_fetcher.get_organism_refname(self.name)
    if not os.path.exists(self.sampledir):
      self.logger.error("Invalid file path to infolder, {}".format(self.sampledir))
      sys.exit()
    try:
      with open("{}".format(infile), 'r') as insample:
        for line in insample:
          #Ignore commented fields
          if not line[0] == '#':
            elem_list = line.rstrip().split("\t")
            if not elem_list[1] == 'N/A':
              hypo.append(dict())
              hypo[-1]["CG_ID_sample"] = self.name
              hypo[-1]["identity"] = elem_list[4]
              hypo[-1]["evalue"] = elem_list[5]
              hypo[-1]["bitscore"] = elem_list[6]
              if elem_list[7] < elem_list[8]:
                hypo[-1]["contig_start"] = int(elem_list[7])
                hypo[-1]["contig_end"] = int(elem_list[8])
              else:
                hypo[-1]["contig_start"] = int(elem_list[8])
                hypo[-1]["contig_end"] = int(elem_list[7])

              hypo[-1]["subject_length"] = int(elem_list[11])
         
              # Split elem 3 in loci (name) and allele (number)
              partials = re.search('(.+)_(\d+){1,3}(?:_(\w+))*', elem_list[3]) 
              hypo[-1]["loci"] = partials.group(1)
              hypo[-1]["allele"] = int(partials.group(2))
              hypo[-1]["span"] = float(hypo[-1]["subject_length"])/self.get_locilength('Seq_types', organism, partials.group(0))

              # split elem 2 into contig node_NO, length, cov
              nodeinfo = elem_list[2].split('_')
              hypo[-1]["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
              hypo[-1]["contig_length"] = int(nodeinfo[3])
              hypo[-1]["contig_coverage"] = nodeinfo[5]
    except Exception as e:
      self.logger.error("{}".format(str(e)))

    #Reduction to top hit
    ind = 0
    while ind < len(hypo)-1:
      targ = ind+1
      while targ < len(hypo):
        ignore = False
        #Rightmost is worse
        if float(hypo[ind]["identity"])*(1-abs(1-hypo[ind]["span"])) >= float(hypo[targ]["identity"])*(1-abs(1-hypo[targ]["span"])):
          #self.logger.warn("Removing {}:{} with span {} and id {}".format(hypo[targ]['loci'],hypo[targ]["allele"], hypo[targ]['span'],hypo[targ]['identity']))
          del hypo[targ]
          ignore = True
        #Leftmost is worse
        elif float(hypo[ind]["identity"])*(1-abs(1-hypo[ind]['span'])) < float(hypo[targ]["identity"])*(1-abs(1-hypo[targ]['span'])):
          #self.logger.warn("Removing {}:{} with span {} and id {}".format(hypo[ind]['loci'],hypo[ind]["allele"], hypo[ind]['span'],hypo[ind]['identity']))
          del hypo[ind]
          targ = ind +1
          ignore = True
        if not ignore:
          targ += 1
        else:
          pass
      ind += 1
    for hit in hypo:
      self.logger.info("Kept {}:{} with span {} and id {}".format(hit['loci'],hit["allele"], hit['span'],hit['identity']))
      self.db_pusher.add_rec(hit, 'Seq_types')

  def scrape_alignment(self):
    """Scrapes a single alignment result"""
    ins_list = list()
    cov_dict = dict()
    align_dict = dict()
    align_dict["reference_genome"] = self.lims_fetcher.data['reference']

    #Reading
    q_list = glob.glob("{}/alignment/*.stats.*".format(self.sampledir))
    map_rate = 0.0
    median_ins = 0
    ref_len = 0.0
    tot_reads = 0
    tot_map = 0
    duprate = 0.0
    for file in q_list:
      with open(file, 'r') as fh:
       type = file.split('.')[-1]
       for line in fh.readlines():
         lsplit = line.rstrip().split('\t')
         if type == 'raw':
           tot_reads = int(lsplit[0])
         elif type == 'ins':
           if len(lsplit) == 18 and lsplit[7] in ['FF','FR']:
             try:
               median_ins = int(lsplit[0])
             except Exception as e:
               pass
         elif type == 'cov':
           cov_dict[lsplit[1]] = int(lsplit[2])
         elif type == 'ref':
           if lsplit[0] != '*' and len(lsplit) >= 2:
             ref_len = ref_len + int(lsplit[1])
         elif type == 'dup':
           if lsplit[0] == 'Unknown Library':
             try:
               duprate = float(lsplit[8]) 
             except Exception as e:
               duprate = -1.0
         elif type == 'map':
           dsplit = line.rstrip().split(' ')
           if len(dsplit)>= 5 and dsplit[4] == 'total':
             tot_map = int(dsplit[0])
           elif len(dsplit)>=4 and dsplit[3] == 'mapped':
             if tot_map > 0:
               map_rate = int(dsplit[0])/float(tot_map)

    #Mangling
    sumz, plus10, plus30, plus50, plus100, total = 0, 0, 0, 0, 0, 0
    for k, v in cov_dict.items():
      sumz += int(k)*v
      total += v
      if int(k) > 10:
        plus10 += v
      if int(k) > 30:
        plus30 += v
      if int(k) > 50:
        plus50 += v
      if int(k) > 100:
        plus100 += v
    if total > 0:
      align_dict['coverage_10x'] = plus10/float(ref_len)
      align_dict['coverage_30x'] = plus30/float(ref_len)
      align_dict['coverage_50x'] = plus50/float(ref_len)
      align_dict['coverage_100x'] = plus100/float(ref_len)
    else:
      align_dict['coverage_10x'] = 0.0
      align_dict['coverage_30x'] = 0.0
      align_dict['coverage_50x'] = 0.0
      align_dict['coverage_100x'] = 0.0
 
    align_dict['mapped_rate'] = map_rate
    align_dict['insert_size'] = median_ins
    if ref_len > 0:
      align_dict['duplication_rate'] = duprate
      align_dict['average_coverage'] = sumz/float(ref_len)
    else:
      align_dict['duplication_rate'] = 0.0
      align_dict['average_coverage'] = 0.0
    align_dict['total_reads'] = tot_reads
    self.db_pusher.upd_rec({'CG_ID_sample' : self.name}, 'Samples', align_dict)
