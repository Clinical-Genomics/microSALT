"""This initial script creates sbatch jobs
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import re
import time
import yaml

import pdb # debug

with open("{}/config.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
  config = yaml.load(conf)

## This needs to be revised
@click.command()
@click.argument('indir')

def main(indir):
  indir = os.path.abspath(indir)
  now = time.strftime("%Y.%m.%d_%H.%M.%S")
  name = os.path.basename(os.path.normpath(indir))
  outdir = "{}/{}_{}".format(config["folders"]["results"],name, now)
  if not os.path.exists(outdir):
    os.makedirs(outdir)
  batchfile = open("{}/runfile.sbatch".format(outdir), "w+")
  #SBatch header
  batchfile.write("#!/bin/bash -l\n\n")
  batchfile.write("#SBATCH -A {}\n".format(config["slurm_header"]["project"]))
  batchfile.write("#SBATCH -p {}\n".format(config["slurm_header"]["type"]))
  batchfile.write("#SBATCH -n {}\n".format(config["slurm_header"]["threads"]))
  batchfile.write("#SBATCH -t {}\n".format(config["slurm_header"]["time"]))
  batchfile.write("#SBATCH -J {}_{}\n".format(config["slurm_header"]["job_name"], now))
  batchfile.write("#SBATCH --qos {}\n\n".format(config["slurm_header"]["qos"]))
  batchfile.close()

  #SPAdes job
  #Loose file format
  fileformat = re.compile('(\d{1}_\d{6}_\w{9}_\w{10,12}_\w{8,12}_)(\d{1})(.fastq.gz)')
  files = os.listdir(indir)
  verified_files = list()
  while len(files) > 0:
    file_parts = fileformat.match( files.pop(0) )
    #If file meets standard format, find pair
    if file_parts:
      if file_parts[2] == '1':
        pairno = '2'
      elif file_parts[2] == '2':
        pairno = '1'
      else:
        print("No pair found. Add an error message.") 

      pairname = "{}{}{}".format(file_parts[1], pairno, file_parts[3])
      if pairname in files:
        files.pop( files.index(pairname) )
        verified_files.append(file_parts[0])
        verified_files.append(pairname)
      else:
        print("No pair found. Add an error message.")

  batchfile = open("{}/runfile.sbatch".format(outdir), "a+")
  #memory is actually 128 per node regardless of cores.
  batchfile.write("spades.py --threads {} --memory {} -o {}/assembly"\
  .format(config["slurm_header"]["threads"], 8*int(config["slurm_header"]["threads"]), outdir))

  pairno = 0
  for file in verified_files:
    #If index is even, raise pair no
    if (verified_files.index(file) & 1) == 0:
      pairno = pairno + 1
    #Set read no
    readno = verified_files.index(file) % 2 + 1
    batchfile.write(" --pe{}-{} {}/{}".format(pairno, readno, indir, file))
  batchfile.write("\n\n")
  batchfile.close()

  #BLAST job
  #index database
  batchfile = open("{}/runfile.sbatch".format(outdir), "a+")
  #Organism needs to be fetched in the future 
  batchfile.write("cd {} && makeblastdb -in {}/{}.xmfa -dbtype nucl -parse_seqids -out {}\n".format(config["folders"]["references"], config["folders"]["references"], config["organism"], config["organism"]))
  #create run
  blast_format = "7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send"
  batchfile.write("blastn -db {}/{} -query {}/assembly/contigs.fasta -out {}/loci_query_tab.txt -num_threads {} -max_target_seqs 1 -outfmt {}\n\n".format(config["folders"]["references"], config["organism"], outdir, outdir, config["slurm_header"]["threads"], blast_format))

	
if __name__ == '__main__':
    main()
