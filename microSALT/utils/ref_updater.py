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

class Ref_Updater():

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

  def fetch_pubmlst(self):
    """ Updates reference for data that is stored on pubMLST """
    # Example request URI: http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1/profiles_csv
    seqdef_url = dict()
    databases = "http://rest.pubmlst.org/db"
    req = urllib.request.Request(databases)
    with urllib.request.urlopen(req) as response:
      query = json.loads(response.read().decode('utf-8'))
    
    # Fetch seqdef locations 
    for name in self.organisms:
      for item in query:
        for subtype in item['databases']:
          if name.replace('_', ' ') in subtype['description'].lower():
            #Seqdef always appear after isolates, so this is fine
            self.updated.append(name.replace('_', ' '))
            seqdef_url[name] = subtype['href']

    for key, val in seqdef_url.items():
      pname = "profile_{}".format(key)
      currver = self.db_access.get_version(pname)

      # Fetch version
      ver_input = "{}/schemes/1/profiles".format(val)
      req = urllib.request.Request(ver_input)
      with urllib.request.urlopen(req) as response:
        query = json.loads(response.read().decode('utf-8'))
      if query['last_updated'] > currver:
        # Download definition files
        output = "{}/{}".format(self.config['folders']['profiles'], key)
        input = "{}/schemes/1/profiles_csv".format(val)
        urllib.request.urlretrieve(input, output)
        # Update database
        self.db_access.upd_rec({'name':pname}, 'Versions', {'version':query['last_updated']})
        self.db_access.reload_profiletable(key)

        loci_input="{}/schemes/1".format(val)
        req = urllib.request.Request(loci_input)
        with urllib.request.urlopen(req) as response:
          query = json.loads(response.read().decode('utf-8'))
        # Clear existing directory and download allele files
        out = "{}/{}".format(self.config['folders']['references'], key)
        shutil.rmtree(out)
        os.makedirs(out)
        # Example request URI: http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/abcZ/alleles_fasta
        for locipath in query['loci']:
          loci = os.path.basename(os.path.normpath(locipath))
          urllib.request.urlretrieve("{}/alleles_fasta".format(locipath), "{}/{}.tfa".format(out, loci))
