"""Scrapes output files for data and adds them to the database
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import glob
import os
import re
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
      target = re.search('(.+)_(\w+)', target).group(1)
      filename="{}/{}/{}.tfa".format(self.config["folders"]["references"], reference, target)

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
                hypo[-1]["contig_start"] = elem_list[7]
                hypo[-1]["contig_end"] = elem_list[8]
              else:
                hypo[-1]["contig_start"] = elem_list[8]
                hypo[-1]["contig_end"] = elem_list[7]
              hypo[-1]["subject_length"] =  elem_list[11]

              # Split elem 3 in loci (name) and allele (number)
              partials = re.search('(.+)_(\d+){1,2}(?:_(\w+))*', elem_list[3])
              hypo[-1]["reference"] = partials.group(3)

              hypo[-1]["gene"] = partials.group(1)
              if hypo[-1]["gene"] in self.gene2resistance.keys():
                hypo[-1]["resistance"] = self.gene2resistance[hypo[-1]["gene"]]
              else:
                #Due to ass backwards naming conventions in recent resFinder, checking backwards
                trial = hypo[-1]["gene"]
                while len(trial) > 0: 
                  trial = trial[:-1]
                  if trial in self.gene2resistance.keys(): 
                    hypo[-1]['gene'] = trial
                    hypo[-1]["resistance"] = self.gene2resistance[trial]
                    break
                if len(trial) == 0:
                  raise Exception("Unable to associate resistance for {} of {} in {}".format(elem_list[3],hypo[-1]['CG_ID_sample'], os.path.basename(file)))
                self.logger.warn("Alternerative split for  hit {} of {}. Resolved to {}".format(elem_list[3], hypo[-1]['CG_ID_sample'], hypo[-1]['gene']))

              hypo[-1]["instance"] = os.path.basename(file[:-4])
              hypo[-1]["span"] = float(hypo[-1]["subject_length"])/self.get_locilength('Resistances', hypo[-1]["instance"], partials.group(0))

              # split elem 2 into contig node_NO, length, cov
              nodeinfo = elem_list[2].split('_')
              hypo[-1]["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
              hypo[-1]["contig_length"] = nodeinfo[3]
              hypo[-1]["contig_coverage"] = nodeinfo[5]

    self.logger.info("{} candidate resistance hits found".format(len(hypo)))

    #Keep track of this section. Bound to goof up
    #Cleanup of overlapping hits
    ind = 0
    while ind < len(hypo)-1:
      targ = ind+1
      while targ < len(hypo):
        ignore = False
        if hypo[ind]["contig_name"] == hypo[targ]["contig_name"]:      
          #Overlapping or shared gene 
          debString = "{} ({}-{}) [id:{},span:{}] overlaps with (or shares name) {} ({}-{}) [id:{},span:{}].".format(hypo[ind]['gene'],hypo[ind]['contig_start'],hypo[ind]['contig_end'],hypo[ind]['identity'],hypo[ind]['span'],hypo[targ]['gene'],hypo[targ]['contig_start'],hypo[targ]['contig_end'],hypo[targ]['identity'],hypo[targ]['span'])
          if (hypo[ind]["contig_start"] <= hypo[targ]["contig_end"] and hypo[ind]["contig_end"] >= hypo[targ]["contig_start"]) or hypo[ind]['gene'] == hypo[targ]['gene']:
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
          #self.logger.info("Debug: {}".format(debString))
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

    #Exceptions
    conversions["aac(6')-Ib3"]="Aminoglycoside"
    conversions["aac(6')-30-aac(6')-Ib'"]="Aminoglycoside"
    conversions["aac(3)-Ib-aac(6')-Ib'"]="Aminoglycoside"
    conversions["aac(6')-Ib-11"]="Aminoglycoside"
    conversions["aac(6')-Ib-Suzhou"]="Aminoglycoside"
    conversions["aac(6')-Ib-Hangzhou"]="Aminoglycoside"
    conversions["aph(3')-VI"]="Aminoglycoside"
    conversions["aph(3')-VIj"]="Aminoglycoside"
    conversions["ant(3'')-Ia"]="Aminoglycoside"
    conversions["aph(3'')-Ib"]='Aminoglycoside'
    conversions['aph(6)-Id']='Aminoglycoside'
    conversions["aac(6')-30"]='Aminoglycoside' 
    conversions['aacA4']='Aminoglycoside'

    conversions["blaSHV-1b-b"]="Beta-lactam"
    conversions['blaSHV-1b']='Beta-lactam'
    conversions['blaCMY-2b']='Beta-lactam'
    conversions['mdf(A)']='Macrolide'
    qnrS = range(1, 10)
    for num in qnrS:
      conversions['qnrS{}'.format(num)] = 'Quinolone'
      conversions['qepA{}'.format(num)] = 'Quinolone'
      conversions['qnrA{}'.format(num)] = 'Quinolone'
    bla = range(1, 300)
    for num in bla:
      conversions['qnrB{}'.format(num)] = 'Quinolone'
      conversions['blaSHV-{}'.format(num)] = 'Beta-lactam'
      conversions['blaCTX-M-{}'.format(num)] = 'Beta-lactam'
      conversions['blaTEM-{}'.format(num)] = 'Beta-lactam'
      conversions['blaVIM-{}'.format(num)] = 'Beta-lactam'
      conversions['blaNDM-{}'.format(num)] = 'Beta-lactam'
      conversions['blaDHA-{}'.format(num)] = 'Beta-lactam'
      conversions['blaCMY-{}'.format(num)] = 'Beta-lactam'

    oxa = range(1, 700)
    for num in oxa:
      conversions['blaOXA-{}'.format(num)] = 'Beta-lactam'
      conversions['blaOXA{}'.format(num)] = 'Beta-lactam'


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
              if elem_list[7] < elem_list[8]:
                seq_col["contig_start"] = elem_list[7]
                seq_col["contig_end"] = elem_list[8]
              else:
                seq_col["contig_start"] = elem_list[8]
                seq_col["contig_end"] = elem_list[7]

              seq_col["contig_start"] = elem_list[7]
              seq_col["contig_end"] = elem_list[8]
              seq_col["subject_length"] =  elem_list[11]
         
              # Split elem 3 in loci (name) and allele (number)
              partials = re.search('(.+)_(\d+){1,2}(?:_(\w+))*', elem_list[3]) 
              seq_col["loci"] = partials.group(1)
              seq_col["allele"] = int(partials.group(2))
              seq_col["span"] = float(seq_col["subject_length"])/self.get_locilength('Seq_types', organism, partials.group(0))

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
