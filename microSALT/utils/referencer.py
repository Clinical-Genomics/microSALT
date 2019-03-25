"""Compares existing organism references with available and updates as needed
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import glob
import json
import os
import re
import shutil
import subprocess
import urllib.request

from Bio import Entrez
from bs4 import BeautifulSoup
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.store.lims_fetcher import LIMS_Fetcher

class Referencer():

  def __init__(self, config, log, force=False):
    self.config = config
    self.logger = log
    self.db_access = DB_Manipulator(config, log)
    self.updated = list()
    #Fetch names of existing refs
    self.refs = self.db_access.profiles
    organisms = self.refs.keys()
    self.organisms = [*organisms]
    self.lims=LIMS_Fetcher(config, log)
    self.force = force

  def identify_new(self, cg_id, project=False):
   """ Automatically downloads pubMLST & NCBI organisms not already downloaded """
   neworgs = list()
   newrefs = list()
   try:
     if project:
       samplenames = self.lims.samples_in_project(cg_id)
       for cg_sampleid in samplenames:
         self.lims.load_lims_sample_info(cg_sampleid)
         refname = self.lims.get_organism_refname(cg_sampleid)
         if refname not in self.organisms and refname not in neworgs:
           neworgs.append(self.lims.data['organism'])
         if not "{}.fasta".format(self.lims.data['reference']) in os.listdir(self.config['folders']['genomes']) and not self.lims.data['reference'] in newrefs:
           newrefs.append(self.lims.data['reference']) 
       for org in neworgs:
         self.add_pubmlst(org)
         #self.download_external(org)
       for org in newrefs:
         self.download_ncbi(org)
     else:
       self.lims.load_lims_sample_info(cg_id)
       refname = self.lims.get_organism_refname(cg_id)
       if refname not in self.organisms:
         self.add_pubmlst(self.lims.data['organism'])
         #self.download_external(self.lims.data['organism'])
       if not "{}.fasta".format(self.lims.data['reference']) in os.listdir(self.config['folders']['genomes']):
         self.download_ncbi(self.lims.data['reference'])
   except Exception as e:
     raise Exception("Unable to add reference for sample {}. Either pubMLST lacks organism {} or NCBI lacks ref {}".format(cg_id, self.lims.data['organism'], self.lims.data['reference']))
 
  def update_refs(self):
    """Updates all references. Order is important, since no object is updated twice"""
    self.fetch_pubmlst(self.force)
    self.fetch_external(self.force)
    self.fetch_resistances(self.force)

  def index_db(self, full_dir, suffix):
    """Check for indexation, makeblastdb job if not enough of them."""
    files = os.listdir(full_dir)
    sufx_files = glob.glob("{}/*{}".format(full_dir, suffix)) #List of source files
    nin_suff = sum([1 for elem in files if 'nin' in elem]) #Total number of nin files
    #if nin_suff < len(sufx_files):
    for file in sufx_files:
      try:
        #Resistence files
        if '.fsa' in suffix:
          bash_cmd = "makeblastdb -in {}/{} -dbtype nucl -out {}".format(\
          full_dir, os.path.basename(file),  os.path.basename(file[:-4]))
        #MLST locis
        else:
          bash_cmd = "makeblastdb -in {}/{} -dbtype nucl -parse_seqids -out {}".format(\
          full_dir, os.path.basename(file),  os.path.basename(file[:-4]))
        proc = subprocess.Popen(bash_cmd.split(), cwd=full_dir, stdout=subprocess.PIPE)
        output, error = proc.communicate()
      except Exception as e:
        self.logger.error("Unable to index requested target {} in {}".format(file, full_dir))
    self.logger.info("Indexed contents of {}".format(full_dir))

  def fetch_external(self, force=False):
    """ Updates reference for data that IS ONLY LINKED to pubMLST """
    prefix = "https://pubmlst.org"
    query = urllib.request.urlopen("{}/data/".format(prefix))
    soup = BeautifulSoup(query, 'html.parser')
    tr_sub = soup.find_all("tr", class_="td1")

    # Only search every other instance
    iterator = iter(tr_sub)
    unfound = True
    try:
      while unfound:
        entry = iterator.__next__()
        # Gather general info from first object
        sample = entry.get_text().split('\n')
        organ = sample[1].lower().replace(' ', '_')
        # In order to get ecoli #1
        if "escherichia_coli" in organ and "#1" in organ:
          organ = organ[:-2]
        currver = self.db_access.get_version("profile_{}".format(organ))
        profile_no = re.search('\d+', sample[2]).group(0)
        if organ in self.organisms and organ.replace("_", " ") not in self.updated and (profile_no > currver or force):
          # Download definition files
          st_link = prefix + entry.find_all("a")[1]['href']
          output = "{}/{}".format(self.config['folders']['profiles'], organ)
          urllib.request.urlretrieve(st_link, output)
          # Update database
          self.db_access.upd_rec({'name':"profile_{}".format(organ)}, 'Versions', {'version':profile_no})
          self.db_access.reload_profiletable(organ)
          # Gather loci from second object
          entry = iterator.__next__()
          # Clear existing directory and download allele files
          out = "{}/{}".format(self.config['folders']['references'], organ)
          shutil.rmtree(out)
          os.makedirs(out)
          for loci in entry.find_all("a"):
            loci = loci['href']
            lociname = os.path.basename(os.path.normpath(loci))
            input = prefix + loci
            urllib.request.urlretrieve(input, "{}/{}".format(out, lociname))
          # Create new indexes
          self.index_db(out, '.tfa')
        else:
          iterator.__next__()
    except StopIteration:
      pass

  def fetch_resistances(self, force=False):
    cwd = os.getcwd()
    url = "https://bitbucket.org/genomicepidemiology/resfinder_db.git"
    hiddensrc ="{}/.resfinder_db".format(self.config['folders']['resistances'])
    wipeIndex = False

    actual = os.listdir(self.config['folders']['resistances'])
    for file in os.listdir(hiddensrc):
      if file not in actual and ('.fsa' in file or 'notes' in file):
        self.logger.info("resFinder database files corrupted. Syncing...")
        wipeIndex = True

    if not os.path.isdir(hiddensrc):
      self.logger.info("resFinder database not found. Fetching..")
      os.makedirs(hiddensrc)
      cmd = "git clone {} --quiet".format(url)
      process = subprocess.Popen(cmd.split(),cwd=self.config['folders']['resistances'], stdout=subprocess.PIPE)
      output, error = process.communicate()
      os.rename("{}/resfinder_db".format(self.config['folders']['resistances']), hiddensrc)
      wipeIndex = True

    else:
      cmd = "git pull origin master"
      process = subprocess.Popen(cmd.split(),cwd=hiddensrc, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
      output, error = process.communicate()
      if not 'Already up-to-date.' in str(output):
        self.logger.info("Remote resFinder database updated. Syncing...")
        wipeIndex = True
      else:
        self.logger.info("Cached resFinder database identical to remote.")

    #Actual update of resistance folder
    if wipeIndex:
      for file in os.listdir(hiddensrc):
        if(os.path.isfile("{}/{}".format(hiddensrc, file))):
          #Copy fresh
          shutil.copy("{}/{}".format(hiddensrc, file), self.config['folders']['resistances'])

    #Double checks indexation is current.
    reIndex = False
    for file in os.listdir(self.config['folders']['resistances']):
      parts = file.split('.')
      if len(parts) > 1 and parts[1] == 'fsa':
        # Missing index or source modified after index
        if not "{}.nhr".format(parts[0]) in os.listdir(self.config['folders']['resistances']) \
           or os.stat("{}/{}".format(self.config['folders']['resistances'], file)).st_mtime > os.stat("{}/{}.nhr".format(self.config['folders']['resistances'],parts[0])).st_mtime\
           or force:
             reIndex = True
    if reIndex:
      self.index_db(self.config['folders']['resistances'], '.fsa')

  def existing_organisms(self):
    """ Returns list of all organisms currently added """
    return self.organisms

  def download_ncbi(self, reference):
    """ Checks available references, downloads from NCBI if not present """
    DEVNULL = open(os.devnull, 'wb')
    Entrez.email="2@2.com"
    record = Entrez.efetch(db='nucleotide', id=reference, rettype='fasta', retmod='text')
    sequence = record.read()
    output = "{}/{}.fasta".format(self.config['folders']['genomes'], reference)
    with open(output, 'w') as f:
      f.write(sequence)
    bwaindex = "bwa index {}".format(output)
    proc = subprocess.Popen(bwaindex.split(), cwd=self.config['folders']['genomes'], stdout=DEVNULL, stderr=DEVNULL)
    out, err = proc.communicate()
    samindex = "samtools faidx {}".format(output)
    proc = subprocess.Popen(samindex.split(), cwd=self.config['folders']['genomes'], stdout=DEVNULL, stderr=DEVNULL)
    out, err = proc.communicate()
    
    self.logger.info('Downloaded reference {}'.format(reference))


  def add_pubmlst(self, organism):
    """ Checks pubmlst for references of given organism and downloads them """
    #Organism must be in binomial format and only resolve to one hit
    organism = organism.lower().replace('.',' ')
    if organism.replace(' ', '_') in self.organisms and not self.force:
      self.logger.info("Organism {} already stored in microSALT".format(organism))
      return
    db_query = self.query_pubmlst()

    #Doublecheck organism name is correct and unique
    orgparts = organism.split(' ')
    counter = 0.0
    for item in db_query:
      for subtype in item['databases']:
        missingPart = False
        for part in orgparts:
          if not part in subtype['description'].lower():
            missingPart = True
        if not missingPart:
          #Seqdef always appear after isolates, so this is fine
          seqdef_url = subtype['href']
          desc = subtype['description']
          counter += 1.0
          self.logger.info("Located pubMLST hit {} for sample".format(desc))
    if counter > 2.0:
      raise Exception("Organism request resolved to {} organisms. Please be more stringent".format(int(counter/2)))
    elif counter < 1.0:
      #add external
      raise Exception("Unable to find requested organism in pubMLST database")  
    else:
      truename = desc.lower().split(' ')
      truename = "{}_{}".format(truename[0], truename[1])
      self.download_pubmlst(truename, seqdef_url)
      #Update organism list
      self.refs = self.db_access.profiles
      self.logger.info("Created table profile_{}".format(truename))

  def query_pubmlst(self):
    """ Returns a json object containing all organisms available via pubmlst.org """
    # Example request URI: http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1/profiles_csv
    seqdef_url = dict()
    databases = "http://rest.pubmlst.org/db"
    db_req = urllib.request.Request(databases)
    with urllib.request.urlopen(db_req) as response:
      db_query = json.loads(response.read().decode('utf-8'))
    return db_query

  def download_pubmlst(self, organism, subtype_href, force=False):
    """ Downloads ST and loci for a given organism stored on pubMLST if it is more recent. Returns update date """
    organism = organism.lower().replace(' ', '_')

    #Pull ST file
    st_target = "{}/{}".format(self.config['folders']['profiles'], organism)
    input = "{}/schemes/1/profiles_csv".format(subtype_href)
    urllib.request.urlretrieve(input, st_target)
    #Pull locus files
    loci_input="{}/schemes/1".format(subtype_href)
    loci_req = urllib.request.Request(loci_input)
    with urllib.request.urlopen(loci_req) as response:
     loci_query = json.loads(response.read().decode('utf-8'))

    output = "{}/{}".format(self.config['folders']['references'], organism)
    if(os.path.isdir(output)):
      shutil.rmtree(output)
    os.makedirs(output)

    for locipath in loci_query['loci']:
      loci = os.path.basename(os.path.normpath(locipath))
      urllib.request.urlretrieve("{}/alleles_fasta".format(locipath), "{}/{}.tfa".format(output, loci))
    # Create new indexes
    self.index_db(output, '.tfa')

  def external_version(self, organism, subtype_href):
    ver_req = urllib.request.Request("{}/schemes/1/profiles".format(subtype_href))
    with urllib.request.urlopen(ver_req) as response:
        ver_query = json.loads(response.read().decode('utf-8'))
    return ver_query['last_updated']

  def fetch_pubmlst(self,force=False):
    """ Updates reference for data that is stored on pubMLST """
    seqdef_url=dict()
    db_query = self.query_pubmlst()
    
    # Fetch seqdef locations 
    for item in db_query:
      for subtype in item['databases']:
        for name in self.organisms:
          if name.replace('_', ' ') in subtype['description'].lower():
            #Seqdef always appear after isolates, so this is fine
            self.updated.append(name.replace('_', ' '))
            seqdef_url[name] = subtype['href']

    for key, val in seqdef_url.items():
      internal_ver = self.db_access.get_version('profile_{}'.format(key))
      external_ver = self.external_version(key, val)  
      if internal_ver < external_ver:
        self.download_pubmlst(key, val, force)
        self.db_access.upd_rec({'name':'profile_{}'.format(key)}, 'Versions', {'version':internal_ver})
        self.db_access.reload_profiletable(key)
        self.logger.info('pubMLST reference for {} set to version {}'.format(key.replace('_',' ').capitalize(), internal_ver))
