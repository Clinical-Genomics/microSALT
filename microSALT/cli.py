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
  click.echo("Execution finished!")

@root.command()
@click.option('--ref_type', type=click.Choice(['ids', 'paths']), help="Use either CG IDs or sample folder paths", required=True)
@click.option('--input_type', type=click.Choice(['list', 'file']), help="Use either list or static file for input", required=True)
@click.option('--file', help="File containing list of samples to analyse, one per row")
@click.option('--sample', '-s', help="List of samples, one per variable invocation", multiple=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def snp(ctx, ref_type, input_type, file, sample, config, email):
  """Pair-wise SNP distance"""
  allpresent = True
  ctx.obj['config']['snp'] = []
  ctx.obj['config']['regex']['mail_recipient'] = email
  if config != '':
    try:
      with open(os.path.abspath(config), 'r') as conf:
        ctx.obj['config'] = json.load(conf)
    except Exception as e:
      pass

  snplist = []
  if input_type=='file':
    sh = open(os.path.abspath(sample), 'r')
    sf = sh.readlines()
    for line in sf:
      snplist.append(line.rstrip())
  elif input_type=='list':
    snplist = sample

  if ref_type == 'paths':
    for line in snplist:
      if not os.path.isdir("{}/alignment".format(line)):
        click.echo("ERROR: {} does not contain an alignment folder + file".format(line))
        sys.exit(-1)
  else:
    scientist=LIMS_Fetcher(ctx.obj['config'], ctx.obj['log'])
    tlist = []
    for line in snplist:
      try:
        scientist.load_lims_sample_info(line)
      except Exception as e:
        click.echo("Unable to load LIMS sample info for {}".format(line))
        sys.exit(-1) 
      project = scientist.data['CG_ID_project']
      samhits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(line))]
      prohits = [x for x in os.listdir(ctx.obj['config']['folders']['results']) if x.startswith("{}_".format(project))]
      if len(samhits) + len(prohits) > 1:
        click.echo("WARNING: Multiple hits for {}. Selecting latest".format(line))
      elif len(samhits) + len(prohits) == 0:
        allpresent = False
        click.echo("{} does not contain an alignment folder + file. Run analysis for sample".format(line))

      if len(samhits) >= 1:
        tlist.append("{}/{}/alignment".format(ctx.obj['config']['folders']['results'], samhits[-1]))
      elif len(prohits) >= 1:
        tlist.append("{}/{}/{}/alignment".format(ctx.obj['config']['folders']['results'], prohits[-1], line))
    snplist = tlist

  if allpresent:
    overlord = Job_Creator(snplist, ctx.obj['config'], ctx.obj['log'])
    overlord.snp_job() 

@root.group()
@click.pass_context
def type(ctx):
  """Standard QC and typing of project/sample"""
  pass

@type.command()
@click.argument('project_id')
@click.option('--input', help='Full path to project folder',default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--qc_only', help="Only runs QC (alignment stats)", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def project(ctx, project_id, input, dry, config, email, qc_only):
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
      click.echo("Path does not contain project id. Exiting.")
      sys.exit(-1)
  else:
    project_dir = "{}/{}".format(ctx.obj['config']['folders']['seqdata'], project_id)
    if not os.path.isdir(project_dir):
      click.echo("Sequence data folder for {} does not exist.".format(project_id))
      sys.exit(-1)

  click.echo("Checking versions of references..")
  try:
    fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
    fixer.identify_new(project_id,project=True)
    fixer.update_refs()
  except Exception as e:
    click.echo("{}".format(e))
  click.echo("Version check done. Creating sbatch jobs")
  manager = Job_Creator(project_dir, ctx.obj['config'], ctx.obj['log'])
  manager.project_job(qc_only=qc_only)
  done() 

@type.command()
@click.argument('sample_id')
@click.option('--input', help='Full path to sample folder', default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--qc_only', help="Only runs QC (alignment stats)", default=False, is_flag=True)
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def sample(ctx, sample_id, input, dry, config, email, qc_only):
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
    fixer = Referencer(ctx.obj['config'], ctx.obj['log'])
    fixer.identify_new(sample_id,project=False) 
    fixer.update_refs()
    click.echo("Version check done. Creating sbatch job")
    worker = Job_Creator(sample_dir, ctx.obj['config'], ctx.obj['log'])
    worker.project_job(single_sample=True, qc_only=qc_only)
  except Exception as e:
    click.echo("Unable to process sample {} due to '{}'".format(sample_id,e))
  done()

@root.group()
@click.pass_context
def util(ctx):
  """ Utilities for specific purposes """
  pass

@util.group()
@click.pass_context
def finish(ctx):
  """Manually upload analysis and generate results"""
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
      click.echo("Path does not contain project id. Exiting.")
      sys.exit(-1)
  else:
    prohits = ctx.obj['config']['folders']['results'].startswith("{}_".format(project_id))
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
  codemonkey.report()
  done()

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

@util.command()
@click.argument('project_name')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--format', default='html', type=click.Choice(['html', 'csv']))
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
