"""Compares existing organism references with available and updates as needed
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import json
import os
import re
import shutil
import urllib.request

from bs4 import BeautifulSoup
from microSALT.store.db_manipulator import DB_Manipulator

class Referencer():

  def __init__(self, config, log):
    self.config = config
    self.logger = log
    self.db_access = DB_Manipulator(config, log)
    self.updated = list()
    #Fetch names of existing refs
    self.refs = self.db_access.get_profiles()
    organisms = self.refs.keys()
    self.organisms = [*organisms] 
    #cgmlst
    self.cgmlstset = list()

  def generate_cgmlst(self, ref_path):
    """ Takes an input fastq and filters it to for an appropriate gene set """
    self.parse_fastq(reference)
    self.filter_codon()
    self.filter_id()
    outname = "{}/{}.gst".format(self.config["folders"]["gene_set"], os.path.basename(ref_path)[:-4])
    output = open(outname, "w+")
    for item in cgmlstset:
      output.write("{}\n{}\n".format(item[0], item[1]))
    output.close()
    self.logger.info("cgMLST written to {} from {}".format(outname, ref_path))

  def parse_fastq(self, ref):
    """ Parse fastq file and dump into necessary structure """
    with open(ref, 'r') as infile:
      lastkey = ""
      for line in infile:
        line = line.rstrip()
        if '>' in line:
          lastkey = line
        else:
          self.cgmlstset.append(dict())
          self.cgmlstset[-1] = (lastkey , line)

    self.logger.info("Sequences initially: {} items".format(len(self.cgmlstset)))

  def filter_codon(self):
    """ Filter fastq based on start/stop codons and hit length """
    dropped = 0
    for line in self.cgmlstset:
      lastkey = ""
      codons = ['ATG', 'TAG', 'TGA', 'TAA']
      rv_codons = ['TAC', 'ATC', 'ACT', 'ATT']
      #Length and codon filters
      if len(line[1]) >= 50 and (line[1][0:3] in codons and line[1][-3:] in codons) or (line[1][0:3] in rv_codons and line[1][-3:] in rv_codons):
        self.cgmlstset.append(dict())
        self.cgmlstset[-1] = (lastkey , line)
      else:
        dropped += 1

    self.logger.info("Sequences remaining after length and codon filter: {} items".format(len(self.cgmlstset)-dropped))
 
  def filter_id(self):
    """ Filter fastq based on identity to other hits """
    refined = list()
    while len(self.cgmlstset) > 1:
      a = self.cgmlstset.pop(0)
      with open('a', 'w') as temporary:
        temporary.write("{}\n{}".format(a[0],a[1]))
      with open('b', 'w') as temporary:
        for item in self.cgmlstset:
          temporary.write("{}\n{}\n".format(item[0], item[1]))
      cmd = "water -asequence a -bsequence b -gapopen 10.0 -gapextend 0.5 -outfile out.water -datafile EDNAFULL"
      subprocess.check_output(cmd.split(), stderr=subprocess.DEVNULL)
      out = open("{}/out.water".format(os.getcwd()), 'r')
      content = out.readlines()
      #Read output
      similar = False
      for outl in content:
        if 'Length' in outl:
          length = int(re.search('\d+', outl).group(0))
        elif 'Identity' in outl:
          id = float(re.search('\(\ *(\d+.\d+)%\)', outl).group(1))
          #Identity filter when id found (always after length)
          if id >= 90 and length >= 100:
            similar = True
            break
          #Reset
          else:
            length = 0
            id = 0
      if not similar:
        refined.append(dict())
        refined[-1]=(a[0], a[1])

    self.logger.info("Sequences after similarity filter: {} items".format(len(refined)))
    self.cgmlstset = refined

  def update_refs(self):
    """Updates all references. Order is important, since no object is updated twice"""
    self.fetch_pubmlst()
    self.fetch_external()

  def fetch_external(self):
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
        if organ in self.organisms and organ.replace("_", " ") not in self.updated and profile_no > currver:
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
        else:
          iterator.__next__()
    except StopIteration:
      pass

  def existing_organisms(self):
    """ Returns list of all organisms currently added """
    return self.organisms

  def add_pubmlst(self, organism):
    """ Checks pubmlst for references of given organism and downloads them """
    #Organism must be in binomial format and only resolve to one hit
    organism = organism.lower().replace(' ', '_')
    if organism in self.organisms:
      self.logger.info("Organism {} already stored in microSALT".format(organism))
      return
    db_query = self.query_pubmlst()

    #Doublecheck organism name is correct and unique
    counter = 0.0 
    for item in db_query:
      for subtype in item['databases']:
        if organism.replace('_', ' ') in subtype['description'].lower():
          #Seqdef always appear after isolates, so this is fine
          seqdef_url = subtype['href']
          counter += 1.0
    if counter > 2:
      self.logger.info("Organism request resolved to {} organisms. Please be more stringent".format(counter/2))
      return
    elif counter < 2:
      self.logger.info("Unable to find requested organism in pubMLST database")
      return
    self.download_pubmlst(organism, seqdef_url)

  def query_pubmlst(self):
    """ Returns a json object containing all organisms available via pubmlst.org """
    # Example request URI: http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1/profiles_csv
    seqdef_url = dict()
    databases = "http://rest.pubmlst.org/db"
    db_req = urllib.request.Request(databases)
    with urllib.request.urlopen(db_req) as response:
      db_query = json.loads(response.read().decode('utf-8'))
    return db_query

  def download_pubmlst(self, organism, subtype_href):
    """ Downloads ST and loci for a given organism stored on pubMLST if it is more recent. Returns update date """
    organism = organism.lower().replace(' ', '_')

    #Pull version
    ver_req = urllib.request.Request("{}/schemes/1/profiles".format(subtype_href))
    with urllib.request.urlopen(ver_req) as response:
        ver_query = json.loads(response.read().decode('utf-8'))
    currver = self.db_access.get_version("profile_{}".format(organism))
    if ver_query['last_updated'] <= currver:
      return currver

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
    return ver_query['last_updated']

  def fetch_pubmlst(self):
    """ Updates reference for data that is stored on pubMLST """
    seqdef_url=dict()
    db_query = self.query_pubmlst()
    
    # Fetch seqdef locations 
    for name in self.organisms:
      for item in db_query:
        for subtype in item['databases']:
          if name.replace('_', ' ') in subtype['description'].lower():
            #Seqdef always appear after isolates, so this is fine
            self.updated.append(name.replace('_', ' '))
            seqdef_url[name] = subtype['href']

    for key, val in seqdef_url.items():
      ver = self.download_pubmlst(key, val)
      self.db_access.upd_rec({'name':'profile_{}'.format(key)}, 'Versions', {'version':ver})
      self.logger.info('Table profile_{} set to version {}'.format(key, ver))
