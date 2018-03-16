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
        if name.replace('_', ' ') in subtype['description'].lower():
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
    req = urllib.request.Request("{}/schemes/1".format(subtype_href))
    with urllib.request.urlopen(req) as response:
        query = json.loads(response.read().decode('utf-8'))
    currver = self.db_access.get_version("profile_{}".format(organism))
    if query['last_updated'] <= currver:
      return currver

    #Pull ST file
    output = "{}/{}".format(self.config['folders']['profiles'], organism)
    input = "{}/schemes/1/profiles_csv".format(subtype_href)
    urllib.request.urlretrieve(input, output)
    #Pull locus files
    shutil.rmtree(output)
    os.makedirs(output)
    for locipath in query['loci']:
          loci = os.path.basename(os.path.normpath(locipath))
          urllib.request.urlretrieve("{}/alleles_fasta".format(locipath), "{}/{}.tfa".format(output, loci))

    profilename = "profile_{}".format(organism)
    if self.db_access.get_version(profilename) == "0":
      #Add database entry
      self.db_access.add_rec('Versions', {'name':profilename,'version':query['last_updated']})
      self.logger.info("Added {} with version {} to microSALT".format(organism, update))
    else:
      self.db_access.upd_rec({'name':profilename}, 'Versions', {'version':query['last_updated']})
      self.logger.info('Updated {} to version {}'.format(profilename, query['last_updated']))
    self.db_access.reload_profiletable(key)
    return query['last_updated']

  def fetch_pubmlst(self):
    """ Updates reference for data that is stored on pubMLST """
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
      self.download_pubmlst(key, val)
