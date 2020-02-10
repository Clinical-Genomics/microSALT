"""This is the main entry point of microSALT.
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import json
import os
import re
import subprocess
import sys
import yaml

from pkg_resources import iter_entry_points
from microSALT import __version__, config, logger, wd
from microSALT.utils.scraper import Scraper
from microSALT.utils.job_creator import Job_Creator
from microSALT.utils.reporter import Reporter
from microSALT.utils.referencer import Referencer
from microSALT.store.lims_fetcher import LIMS_Fetcher

if config == '':
  click.echo("ERROR - No properly set-up config under neither envvar MICROSALT_CONFIG nor ~/.microSALT/config.json. Exiting.")
  ctx.abort()

def set_cli_config(config):
  if config != '':
    if os.path.exists(config):
      try:
        t = ctx.obj['config']
        with open(os.path.abspath(config), 'r') as conf:
          ctx.obj['config'] = json.load(conf)
        ctx.obj['config']['folders']['expec'] = t['folders']['expec']
        ctx.obj['config']['folders']['adapters'] = t['folders']['adapters']
        ctx.obj['config']['config_path'] = os.path.abspath(config)
      except Exception as e:
        pass

def done():
  click.echo("INFO - Execution finished!")

def validate_param(pfile):
  pass

def pad_param(pfile):
  """Reads the provided parameters json and adds default values as necessary"""

  defaults = {
  "CG_ID_project" : "XXX0000",
  "CG_ID_sample" : "XXX0000Y1",
  "Customer_ID_sample" : "XXX0000Y1",
  "customer_ID" : "cust000",
  "application_tag" : "NONE",
  "date_arrival" : "0001-01-01 00:00:00",
  "date_libprep" : "0001-01-01 00:00:00",
  "date_sequencing" : "0001-01-01 00:00:00",
  "method_libprep" : "Not in LIMS",
  "method_sequencing" : "Not in LIMS",
  "organism" : "Unset",
  "priority" : "standard",
  "reference" : "None"}

  for entry in pfile.items():
    sample_counter = 1
    for k, v in defaults.items():
      if not defaults[k] in entry:
        click.echo("INFO - Parameter {} was not provided. Used default.".format(k)) 
        if k in ['CG_ID_sample','Customer_ID_sample']:
          entry[k] = "XXX0000Y{}".format(sample_counter)
          sample_counter = sample_counter + 1
        else:
          entry[k] = v
      

@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
  """microbial Sequence Analysis and Loci-based Typing (microSALT) pipeline """
  ctx.obj = {}
  ctx.obj['config'] = config
  ctx.obj['log'] = logger

@root.command()
@click.option('--param','-p',default={},help='Json file describing input samples')
@click.option('--input', help='Full path to project folder',default="")
@click.option('--track', help='Run a specific analysis track',default="default", type=click.Choice(['default','typing','qc','cgmlst']))
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--untrimmed', help="Use untrimmed input data", default=False, is_flag=True)
@click.option('--uncareful', help="Avoids running SPAdes in careful mode. Sometimes fix assemblies", default=False, is_flag=True)
@click.pass_context
def analyse(ctx, param, input, track, config, dry, email, skip_update, untrimmed, uncareful):
  """Sequence analysis, typing and resistance identification"""
  #Run section
  pool = []
  trimmed = not untrimmed
  careful = not uncareful
  set_cli_config(config)
  ctx.obj['config']['regex']['mail_recipient'] = email
  ctx.obj['config']['dry'] = dry
  if not os.path.isdir(input):
    click.echo("ERROR - Sequence data folder {} does not exist.".format(input))
    ctx.abort()
  for subfolder in os.listdir(input):
    if os.path.isdir("{}/{}".format(input, subfolder)):
      pool.append(subfolder)

  #Samples section
  validate_param(param)
  param = pad_param(param)
  run_creator = Job_Creator(project_dir, ctx.obj['config'], ctx.obj['log'], param, track=track,trim=trimmed,\
                            careful=careful, pool=pool)

  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("INFO - Checking versions of references..")
  try:
    if not skip_update:
      ext_refs.identify_new(project_id,project=True)
      ext_refs.update_refs()
      click.echo("INFO - Version check done. Creating sbatch jobs")
    else:
      click.echo("INFO - Skipping version check.")
  except Exception as e:
    click.echo("{}".format(e))

  if len(param.items()) > 1:
    run_creator.project_job()
  elif len(param.items()) == 1: 
    run_creator.project_job(single_sample=True)
  else:
    ctx.abort()
 
  done()

@root.group()
@click.pass_context
def utils(ctx):
  """ Utilities for specific purposes """
  pass

@utils.group()
@click.pass_context
def refer(ctx):
  """ Manipulates MLST organisms """
  pass

@utils.command()
@click.option('--param','-p',default={},help='Json file describing input samples')
@click.option('--input', help='Full path to project folder',default="")
@click.option('--track', help='Run a specific analysis track',default="default", type=click.Choice(['default','typing','qc','cgmlst']))
@click.option('--config', help="microSALT config to override default", default="")
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--report', default='default', type=click.Choice(['default', 'typing', 'motif_overview', 'qc', 'json_dump', 'st_update']))
@click.pass_context
def finish(ctx, param, input, track, config, dry, email, skip_update, report):
  """Sequence analysis, typing and resistance identification"""
  #Run section
  pool = []
  set_cli_config(config)
  ctx.obj['config']['regex']['mail_recipient'] = email
  ctx.obj['config']['dry'] = dry
  ctx.obj['config']['rerun'] = True
  if not os.path.isdir(input):
    click.echo("ERROR - Sequence data folder {} does not exist.".format(input))
    ctx.abort()
  for subfolder in os.listdir(input):
    if os.path.isdir("{}/{}".format(input, subfolder)):
      pool.append(subfolder)

  #Samples section
  validate_param(param)
  param = pad_param(param)
  res_scraper = Scraper(project_dir, ctx.obj['config'], ctx.obj['log'], param)

  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("INFO - Checking versions of references..")
  try:
    if not skip_update:
      ext_refs.identify_new(project_id,project=True)
      ext_refs.update_refs()
      click.echo("INFO - Version check done. Creating sbatch jobs")
    else:
      click.echo("INFO - Skipping version check.")
  except Exception as e:
    click.echo("{}".format(e))

  if len(param.items()) >= 1:
    #res_scraper.scrape_project()
    for subfolder in pool:
      res_scraper = Scraper("{}/{}".format(input, subfolder), ctx.obj['config'], ctx.obj['log'])
      res_scraper.scrape_sample()
  else:
    ctx.abort()
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], param, output=input, collection=True)
  codemonkey.report(report)
  done()

@refer.command()
@click.argument('organism')
@click.option('--force', help="Redownloads existing organism", default=False, is_flag=True)
@click.pass_context
def add(ctx, organism, force):
  """ Adds a new internal organism from pubMLST """
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  try:
    referee.add_pubmlst(organism)
  except Exception as e:
    click.echo(e.args[0])
    ctx.abort()
  click.echo("INFO - Checking versions of all references..")
  referee = Referencer(ctx.obj['config'], ctx.obj['log'],force=force)
  referee.update_refs()

@refer.command()
@click.pass_context
def list(ctx):
  """ Lists all stored organisms """
  refe = Referencer(ctx.obj['config'], ctx.obj['log'])
  click.echo("INFO - Currently stored organisms:")
  for org in sorted(refe.existing_organisms()):
    click.echo(org.replace("_"," ").capitalize())

@utils.command()
@click.argument('project_name')
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.option('--type', default='default', type=click.Choice(['default', 'typing', 'motif_overview', 'qc', 'json_dump', 'st_update']))
@click.option('--output',help='Full path to output folder',default="")
@click.option('--collection',default=False, is_flag=True)
@click.pass_context
def report(ctx, project_name, email, type, output, collection):
  """Re-generates report for a project"""
  ctx.obj['config']['regex']['mail_recipient'] = email
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'], project_name, output, collection=collection)
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

@utils.command()
@click.option('--dry', help="Builds instance without posting to SLURM", default=False, is_flag=True)
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def autobatch(ctx, dry, skip_update, email):
  """Analyses all currently unanalysed projects in the seqdata folder"""
  #Trace exists by some samples having pubMLST_ST filled in. Make trace function later
  ctx.obj['config']['regex']['mail_recipient'] = email
  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  if not skip_update:
    ext_refs.update_refs()
    ext_refs.resync()

  process = subprocess.Popen('squeue --format="%50j" -h -r'.split(), stdout=subprocess.PIPE)
  run_jobs, error = process.communicate()
  run_jobs = run_jobs.splitlines()
  run_jobs = [ re.search("\"(.+)\"", str(jobname)).group(1).replace(' ', '') for jobname in run_jobs];
  old_jobs = os.listdir(ctx.obj['config']['folders']['results'])

  for foldah in os.listdir(ctx.obj['config']['folders']['seqdata']):
      #Project name not found in slurm list 
      if len([job for job in run_jobs if foldah in job]) == 0:
        #Project name not found in results directories
        if len([job for job in old_jobs if foldah in job]) == 0:
          if dry:
            click.echo("DRY - microSALT analyse {} --skip_update".format(foldah))
          else:
            process = subprocess.Popen("microSALT analyse project {} --skip_update".format(foldah).split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
        elif dry:
          click.echo("INFO - Skipping {} due to existing analysis in results folder".format(foldah))
      elif dry:
        click.echo("INFO - Skipping {} due to concurrent SLURM run".format(foldah))

@utils.command()
@click.option('--input', help='Full path to project folder',default="")
@click.pass_context
def paramgenerate(ctx):
  """Generates the necessary parameter file for analysis"""
  input_folder = os.path.basename(input)

  defaults = {
  "CG_ID_project" : "XXX0000",
  "CG_ID_sample" : "XXX0000Y1",
  "Customer_ID_sample" : "XXX0000Y1",
  "customer_ID" : "cust000",
  "application_tag" : "NONE",
  "date_arrival" : "0001-01-01 00:00:00",
  "date_libprep" : "0001-01-01 00:00:00",
  "date_sequencing" : "0001-01-01 00:00:00",
  "method_libprep" : "Not in LIMS",
  "method_sequencing" : "Not in LIMS",
  "organism" : "Unset",
  "priority" : "standard",
  "reference" : "None"}

  pool = []
  if not os.path.isdir(input):
    click.echo("ERROR - Sequence data folder {} does not exist.".format(input_folder))
    ctx.abort()
  for subfolder in os.listdir(input):
    if os.path.isdir("{}/{}".format(input, subfolder)):
      pool.append(defaults)
      pool[-1]['CG_ID_project'] = input_folder
      pool[-1]['CG_ID_sample'] = subfolder

  output = open("{}/{}".format(os.getcwd(), input_folder), 'w')
  output.write(pool)
  output.close()
  click.echo("INFO - Created {}.json".format(input_folder))
  done()

@resync.command()
@click.option('--type', default='list', type=click.Choice(['report', 'list']), help="Output format")
@click.option('--customer', default='all', help="Customer id filter")
@click.option('--skip_update', default=False, help="Skips downloading of references", is_flag=True)
@click.option('--email', default=config['regex']['mail_recipient'], help='Forced e-mail recipient')
@click.pass_context
def review(ctx, type, customer, skip_update, email):
  """Generates information about novel ST"""
  #Trace exists by some samples having pubMLST_ST filled in. Make trace function later
  ctx.obj['config']['regex']['mail_recipient'] = email
  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  if not skip_update:
    ext_refs.update_refs()
    ext_refs.resync()
  click.echo("INFO - Version check done. Generating output")
  if type=='report':
    codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
    codemonkey.report(type='st_update', customer=customer)
  elif type=='list':
    ext_refs.resync(type=type)
  done()

@resync.command()
@click.argument('sample_name')
@click.option('--force', default=False, is_flag=True, help="Resolves sample without checking for pubMLST match")
@click.pass_context
def overwrite(ctx,sample_name, force):
  """Flags sample as resolved"""
  ext_refs = Referencer(ctx.obj['config'], ctx.obj['log'])
  ext_refs.resync(type='overwrite', sample=sample_name, ignore=force)
  done()
