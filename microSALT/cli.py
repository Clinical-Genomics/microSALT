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
  click.echo("ERROR: No properly set-up config under neither envvar MICROSALT_CONFIG nor ~/.microSALT/config.json. Exiting.")
  sys.exit(-1)

def done():
  click.echo("Execution finished!")

@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
  """microbial Sequence Analysis and Loci-based Typing (microSALT) pipeline """
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

@root.group()
@click.pass_context
def analyse(ctx):
  """Basic sequence and resistance typing"""
  pass

@root.group()
@click.pass_context
def utils(ctx):
  """ Utilities for specific purposes """
  pass

@utils.group()
@click.pass_context
def finish(ctx):
  """Reupload typing analysis and generate results"""
  pass

@utils.group()
@click.pass_context
def refer(ctx):
  """ Manipulates MLST organisms """
  pass

@analyse.command()
@click.argument('project_id')
@click.option('--input', help='Full path to project folder',default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--qc_only', help="Only runs QC (alignment stats)", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--untrimmed', help="Use untrimmed input data", default=False, is_flag=True)
@click.option('--uncareful', help="Avoids running SPAdes in careful mode. Sometimes fix assemblies", default=False, is_flag=True)
@click.pass_context
def project(ctx, project_id, input, dry, config, email, qc_only, untrimmed, skip_update, uncareful):
  """Analyse a whole project"""
  ctx.obj['config']['regex']['mail_recipient'] = email
  trimmed = not untrimmed
  careful = not uncareful
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
      click.echo("Path does not contain project id. Exiting.")
      sys.exit(-1)
  else:
    project_dir = "{}/{}".format(ctx.obj['config']['folders']['seqdata'], project_id)
    if not os.path.isdir(project_dir):
      click.echo("Sequence data folder for {} does not exist.".format(project_id))
      sys.exit(-1)

  click.echo("Checking versions of references..")
  try:
    if not skip_update:
      fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
      fixer.identify_new(project_id,project=True)
      fixer.update_refs()
      print("Version check done. Creating sbatch jobs")
    else:
      print("Skipping version check.")
  except Exception as e:
    print("{}".format(e))

  manager = Job_Creator(project_dir, ctx.obj['config'], ctx.obj['log'],trim=trimmed,qc_only=qc_only, careful=careful)
  manager.project_job()
  done() 

@analyse.command()
@click.argument('sample_id')
@click.option('--input', help='Full path to sample folder', default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--qc_only', help="Only runs QC (alignment stats)", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--untrimmed', help="Use untrimmed input data", default=False, is_flag=True)
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--uncareful', help="Avoids running SPAdes in careful mode. Sometimes fix assemblies", default=False, is_flag=True)
@click.pass_context
def sample(ctx, sample_id, input, dry, config, email, qc_only, untrimmed, skip_update, uncareful):
  """Analyse a single sample"""
  ctx.obj['config']['regex']['mail_recipient'] = email
  trimmed = not untrimmed
  careful = not uncareful
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
    click.echo("Unable to load LIMS sample info.")
    sys.exit(-1)

  if input != "":
    sample_dir = os.path.abspath(input)
    if not sample_id in sample_dir:
      click.echo("Path does not contain sample id. Exiting.")
      sys.exit(-1)
  else:
    sample_dir = "{}/{}/{}".format(ctx.obj['config']['folders']['seqdata'], scientist.data['CG_ID_project'] ,sample_id)
    if not os.path.isdir(sample_dir):
      click.echo("Sequence data folder for {} does not exist.".format(sample_id))
      sys.exit(-1)

  click.echo("Checking versions of references..")
  try:
    if not skip_update:
      fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
      fixer.identify_new(sample_id,project=False) 
      fixer.update_refs()
      print("Version check done. Creating sbatch job")
    else:
      print("Skipping version check.")
    worker = Job_Creator(sample_dir, ctx.obj['config'], ctx.obj['log'], trim=trimmed,qc_only=qc_only, careful=careful)
    worker.project_job(single_sample=True)
  except Exception as e:
    click.echo("Unable to process sample {} due to '{}'".format(sample_id,e))
  done()

@finish.command()
@click.argument('sample_id')
@click.option('--rerun', is_flag=True, default=False, help='Overwrite existing data')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--input', help='Full path to result sample folder', default="")
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--report', default='default', type=click.Choice(['default', 'qc']))
@click.pass_context
def sample(ctx, sample_id, rerun, email, input, config, report):
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
      click.echo("Path does not contain sample id. Exiting.")
      sys.exit(-1)
  else:
    prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(sample_id))]
    if len(prohits) > 1:
      click.echo("Multiple instances of that analysis exists. Specify full path using --input")
      sys.exit(-1)
    elif len(prohits) <1:
      click.echo("No analysis folder prefixed by {} found.".format(project_id))
      sys.exit(-1)
    else:
      sample_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], prohits[-1])

  scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
  try:
    scientist.load_lims_sample_info(sample_id)
  except Exception as e:
    click.echo("Unable to load LIMS sample info.")
    sys.exit(-1)

  garbageman = Scraper(sample_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_sample()

  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], scientist.data['CG_ID_project'])
  codemonkey.report(report)
  done()

@finish.command()
@click.argument('project_id')
@click.option('--rerun', is_flag=True, default=False, help='Overwrite existing data')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--input', help='Full path to result project folder', default="")
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--report', default='default', type=click.Choice(['default', 'qc']))
@click.pass_context
def project(ctx, project_id, rerun, email, input, config, report):
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
      click.echo("Path does not contain project id. Exiting.")
      sys.exit(-1)
  else:
    prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(project_id))]
    if len(prohits) > 1:
      click.echo("Multiple instances of that analysis exists. Specify full path using --input")
      sys.exit(-1)
    elif len(prohits) <1:
      click.echo("No analysis folder prefixed by {} found.".format(project_id))
      sys.exit(-1)
    else:
      project_dir = "{}/{}".format(ctx.obj['config']['folders']['results'], prohits[-1])

  garbageman = Scraper(project_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_project()
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], project_id)
  codemonkey.report(report)
  done()

@refer.command()
@click.argument('organism')
@click.pass_context
@click.option('--force', help="Redownloads existing organism", default=False, is_flag=True)
def add(ctx, organism, force):
  """ Adds a new internal organism from pubMLST """
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  try:
    referee.add_pubmlst(organism)
  except Exception as e:
    click.echo(e.args[0])
    sys.exit(-1)
  click.echo("Checking versions of all references..")
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  referee.update_refs()

@refer.command()
@click.pass_context
def list(ctx):
  """ Lists all stored organisms """
  refe = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("Currently stored organisms:")
  for org in sorted(refe.existing_organisms()):
    click.echo(org.replace("_"," ").capitalize())

@utils.command()
@click.argument('project_name')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--type', default='default', type=click.Choice(['default', 'typing', 'resistance_overview', 'qc', 'json_dump', 'st_update']))
@click.pass_context
def report(ctx, project_name, email, type):
  """Re-generates report for a project"""
  ctx.obj['config']['regex']['mail_recipient'] = email
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], project_name)
  codemonkey.report(type)
  done()

@utils.command()
@click.pass_context
def view(ctx):
  """Starts an interactive webserver for viewing"""
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
  codemonkey.start_web()

@utils.group()
@click.pass_context
def resync(ctx):
  """Updates internal ST with pubMLST equivalent"""

@resync.command()
@click.option('--type', default='html', type=click.Choice(['report', 'list']), help="Output format")
@click.option('--customer', default='all', help="Customer id filter")
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.pass_context
def review(ctx, type, customer, skip_update):
  """Generates information about novel ST"""
  #Trace exists by some samples having pubMLST_ST filled in. Make trace function later
  fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
  if not skip_update:
    fixer.update_refs()
    fixer.resync()
  print("Version check done. Generating output")
  if type=='report':
    codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
    codemonkey.report(type='st_update', customer=customer)
  elif type=='list':
    fixer.resync(type=type)
  done()

@resync.command()
@click.argument('sample_name')
@click.pass_context
def overwrite(ctx,sample_name):
  """Flags sample as resolved"""
  fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
  fixer.resync(type='overwrite', sample=sample_name)
  done()
