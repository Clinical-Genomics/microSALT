"""This is the main entry point of microSALT.
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import logging
import json
import os
import sys
import yaml

from pkg_resources import iter_entry_points
from microSALT import __version__, config
from microSALT.utils.scraper import Scraper
from microSALT.utils.job_creator import Job_Creator
from microSALT.utils.reporter import Reporter
from microSALT.utils.referencer import Referencer
from microSALT.store.lims_fetcher import LIMS_Fetcher

if config == '':
  print("ERROR: No properly set-up config under neither envvar MICROSALT_CONFIG nor ~/.microSALT/config.json. Exiting.")
  sys.exit(-1)

@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
  """ microbial Sequence Analysis and Loci-based Typing (microSALT) pipeline """
  ctx.obj = {}
  ctx.obj['config'] = config
  logger = logging.getLogger('main_logger')
  logger.setLevel(logging.DEBUG)
  fh = logging.FileHandler(os.path.expanduser(config['folders']['log_file']))
  fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
  logger.addHandler(fh)
  ch = logging.StreamHandler()
  ch.setLevel(logging.INFO)
  ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
  logger.addHandler(ch)
  ctx.obj['log'] = logger

def done():
  print("Execution finished!")

@root.group()
@click.pass_context
def start(ctx):
  """Starts analysis of project/sample"""
  pass

@start.command()
@click.argument('project_id')
@click.option('--input', help='Full path to project folder',default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def project(ctx, project_id, input, dry, config, email):
  """Analyze a whole project"""
  ctx.obj['config']['regex']['mail_recipient'] = email
  if config != '':
    try:
      with open(os.path.abspath(config), 'r') as conf:
        ctx.obj['config'] = json.load(conf)
    except Exception as e:
      pass

  ctx.obj['config']['dry'] = dry
  if input != "":
    project_dir = os.path.abspath(input)
    if not project_id in project_dir:
      print("Path does not contain project id. Exiting.")
      sys.exit(-1)
  else:
    project_dir = "{}/{}".format(ctx.obj['config']['folders']['seqdata'], project_id)
    if not os.path.isdir(project_dir):
      print("Sequence data folder for {} does not exist.".format(project_id))
      sys.exit(-1)

  print("Checking versions of references..")
  try:
    fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
    fixer.identify_new(project_id,project=True)
    fixer.update_refs()
  except Exception as e:
    print("{}".format(e))

  print("Version check done. Creating sbatch jobs")
  manager = Job_Creator(project_dir, ctx.obj['config'], ctx.obj['log'])
  manager.project_job()
  done() 

@start.command()
@click.argument('sample_id')
@click.option('--input', help='Full path to sample folder', default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def sample(ctx, sample_id, input, dry, config, email):
  """Analyze a single sample"""
  ctx.obj['config']['regex']['mail_recipient'] = email
  if config != '':
    try:
      with open(os.path.abspath(config), 'r') as conf:
        ctx.obj['config'] = json.load(conf)
    except Exception as e:
      pass

  ctx.obj['config']['dry'] = dry
  scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  try:
    scientist.load_lims_sample_info(sample_id)
  except Exception as e:
    print("Unable to load LIMS sample info.")
    sys.exit(-1)

  if input != "":
    sample_dir = os.path.abspath(input)
    if not sample_id in sample_dir:
      print("Path does not contain sample id. Exiting.")
      sys.exit(-1)
  else:
    sample_dir = "{}/{}/{}".format(ctx.obj['config']['folders']['seqdata'], scientist.data['CG_ID_project'] ,sample_id)
    if not os.path.isdir(sample_dir):
      print("Sequence data folder for {} does not exist.".format(sample_id))
      sys.exit(-1)

  print("Checking versions of references..")
  try:
    fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
    fixer.identify_new(sample_id,project=False) 
    fixer.update_refs()
    print("Version check done. Creating sbatch job")
    worker = Job_Creator(sample_dir, ctx.obj['config'], ctx.obj['log'])
    worker.project_job(single_sample=True)
  except Exception as e:
    print("Unable to process single sample {} due to '{}'".format(sample_id,e))
  done()

@root.group()
@click.pass_context
def finish(ctx):
  """Uploads analysis and generates results"""
  pass

@finish.command()
@click.argument('sample_id')
@click.option('--rerun', is_flag=True, default=False, help='Overwrite existing data')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--input', help='Full path to result sample folder', default="")
@click.option('--config', help="microSALT config to override default", default="")
@click.pass_context
def sample(ctx, sample_id, rerun, email, input, config):
  """Parse results from analysing a single sample"""
  if config != '':
    try:
      with open(os.path.abspath(config), 'r') as conf:
        ctx.obj['config'] = json.load(conf)
    except Exception as e:
      pass

  ctx.obj['config']['rerun'] = rerun
  ctx.obj['config']['regex']['mail_recipient'] = email
  
  if input != "":
    sample_dir = os.path.abspath(input)
    if not sample_id in sample_dir:
      print("Path does not contain sample id. Exiting.")
      sys.exit(-1)
  else:
    hits = 0
    for i in os.listdir(ctx.obj['config']['folders']['results']):
      if '{}_'.format(sample_id) in i:
        hits = hits+1
        fname = i
    if hits > 1: #Doublechecks only 1 analysis exists
      print("Multiple instances of that analysis exists. Specify full path using --input")
      sys.exit(-1)
    elif hits < 1:
      print("No analysis folder prefixed by {} found.".format(sample_id))
      sys.exit(-1)
    else:
      sample_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], fname)

  scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  try:
    scientist.load_lims_sample_info(sample_id)
  except Exception as e:
    print("Unable to load LIMS sample info.")
    sys.exit(-1)

  garbageman = Scraper(sample_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_sample()

  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], scientist.data['CG_ID_project'])
  codemonkey.report()
  done()

@finish.command()
@click.argument('project_id')
@click.option('--rerun', is_flag=True, default=False, help='Overwrite existing data')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--input', help='Full path to result project folder', default="")
@click.option('--config', help="microSALT config to override default", default="")
@click.pass_context
def project(ctx, project_id, rerun, email, input, config):
  """Parse results from analysing a single project"""
  if config != '':
    try:
      with open(os.path.abspath(config), 'r') as conf:
        ctx.obj['config'] = json.load(conf)
    except Exception as e:
      pass

  ctx.obj['config']['rerun'] = rerun
  ctx.obj['config']['regex']['mail_recipient'] = email
 
  if input != "":
    project_dir = os.path.abspath(input)
    if not project_id in project_dir:
      print("Path does not contain project id. Exiting.")
      sys.exit(-1)
  else:
    hits = 0
    for i in os.listdir(ctx.obj['config']['folders']['results']):
      if '{}_'.format(project_id) in i:
        hits = hits+1
        fname = i
    if hits > 1: #Doublechecks only 1 analysis exists
      print("Multiple instances of that analysis exists. Specify full path using --input")
      sys.exit(-1)
    elif hits < 1:
      print("No analysis folder prefixed by {} found.".format(project_id))
      sys.exit(-1)
    else:
      project_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], fname)
  
  garbageman = Scraper(project_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_project()
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], project_id)
  codemonkey.report()
  done()

@root.group()
@click.pass_context
def util(ctx):
  """ Utilities for specific purposes """
  pass

@util.group()
@click.pass_context
def refer(ctx):
  """ Manipulates MLST organisms """
  pass

@refer.command()
@click.argument('organism')
@click.pass_context
@click.option('--force', help="Redownloads existing organism", default=False, is_flag=True)
def add(ctx, organism, force):
  """ Adds a new internal organism from pubMLST """
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  referee.add_pubmlst(organism)
  print("Checking versions of all references..")
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  referee.update_refs()

@refer.command()
@click.pass_context
def list(ctx):
  """ Lists all stored organisms """
  refe = Referencer(ctx.obj['config'], ctx.obj['log'])
  print("Currently stored organisms:")
  for org in sorted(refe.existing_organisms()):
    print(org.replace("_"," ").capitalize())

@util.command()
@click.argument('project_name')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--format', default='html', type=click.Choice(['html', 'csv', 'json', 'st']))
@click.pass_context
def report(ctx, project_name, email, format):
  """Re-generates report for a project"""
  ctx.obj['config']['regex']['mail_recipient'] = email
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], project_name)
  codemonkey.report(format)
  done()

@util.command()
@click.pass_context
def view(ctx):
  """Starts an interactive webserver for viewing"""
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
  codemonkey.start_web()

@util.group()
@click.pass_context
def resync(ctx):
  """Updates internal ST with pubMLST equivalent"""

@resync.command()
@click.option('--format', default='html', type=click.Choice(['html', 'list']))
@click.option('--customer', default='all')
@click.pass_context
def review(ctx, format, customer):
  """Generates information about novel ST"""
  fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
  fixer.identify_new(sample_id,project=False)
  fixer.update_refs()
  fixer.resync()
  print("Version check done. Generating output")
  if format=='html':
    codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
    codemonkey.report(type='st', customer=customer)
  elif format=='list':
    fixer.resync(type=format)
  done()

@resync.command()
@click.argument('sample_name')
@click.pass_context
def overwrite(ctx,sample_name):
  """All ST with pubMLST equivalent will be marked as resolved"""
  fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
  fixer.resync(type='overwrite')
  done()
