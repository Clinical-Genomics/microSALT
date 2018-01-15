"""Compares existing organism references with available and updates as needed
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import json
import os
import shutil
import urllib.request

from microSALT.store.db_manipulator import DB_Manipulator

class Ref_Updater():

  def __init__(self, config, log):
    self.config = config
    self.logger = log
    self.db_access = DB_Manipulator(config, log)

    self.refs = self.db_access.get_profiles()
    #Fetch names of existing refs
    organisms = self.refs.keys()
    self.organisms = [*organisms]
  
  def update_refs(self):
    self.fetch_pubmlst()

  def fetch_enterobase(self):
    pass

  def fetch_pasteur(self):
    pass

  def fetch_pubmlst(self):
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
            seqdef_url[name] = subtype['href']

    for key, val in seqdef_url.items():
      pname = "profile_{}".format(key)
      currver = self.db_access.get_version(pname)

      # Fetch version
      ver_input = "{}/schemes/1/profiles".format(val)
      req = urllib.request.Request(ver_input)
      with urllib.request.urlopen(req) as response:
        query = json.loads(response.read().decode('utf-8'))
      if query['last_updated'] >= currver:
        self.db_access.upd_rec({'name':pname}, 'Versions', {'version':query['last_updated']})
        # Download definition files
        output = "{}/{}".format(self.config['folders']['profiles'], key)
        input = "{}/schemes/1/profiles_csv".format(val)
        urllib.request.urlretrieve(input, output)
        # Update actual profiletable
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
