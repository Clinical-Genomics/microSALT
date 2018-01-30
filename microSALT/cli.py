"""This is the main entry point of microSALT.
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import json
import os
import sys
import yaml
import logging

from pkg_resources import iter_entry_points
from microsalt import __version__, app
from microsalt.utils.scraper import Scraper
from microsalt.utils.job_creator import Job_Creator
from microsalt.utils.ref_updater import Ref_Updater

# Load configuration
defaulto = os.path.join(os.environ['HOME'], '.microSALT/config.json')
if os.path.exists(defaulto):
  try:
    with open(default, 'r') as conf:
      config = json.load(conf)
  except Exception as e:
    print("ERROR: Config under default path ~/.microSALT/config.json improperly formatted. Exiting")
elif 'MICROSALT_CONFIG' in os.environ:
  try:
    envvar = os.environ['MICROSALT_CONFIG']
    with open(envvar, 'r') as conf:
      config = json.load(conf)
  except Exception as e:
    print("ERROR: Config under envvar MICROSALT_CONFIG improperly formatted. Exiting")
else:
  print("ERROR: No properly set-up config under neither ~/.microSALT/config.json nor envvar MICROSALT_CONFIG. Exiting")
  sys.exit(-1)

# Load flask instance
try:
  app.config.update(config['database'])
except Exception as e:
  logger.error("Config.json lacks 'database' section. Unable to call mySQL. Exiting")
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
  fh = logging.FileHandler(os.path.join(os.environ['HOME'], "microSALT.log"))
  fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
  logger.addHandler(fh)
  ch = logging.StreamHandler()
  ch.setLevel(logging.WARNING)
  ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
  logger.addHandler(ch)
  ctx.obj['log'] = logger

def done():
  print("Execution finished!")

@root.group()
@click.pass_context
def create(ctx):
  """Produces sbatch jobs of the given input"""
  pass

@create.command()
@click.argument('project_dir')
@click.pass_context
def project(ctx, project_dir):
  """Create jobs for a project"""
  print("Checking versions of references..")
  fixer = Ref_Updater(ctx.obj['config'], ctx.obj['log'])
  fixer.update_refs()
  print("Version check done. Creating sbatch jobs")
  manager = Job_Creator(project_dir, ctx.obj['config'], ctx.obj['log'])
  manager.project_job()
  done() 

@create.command()
@click.argument('sample_dir')
@click.pass_context
def sample(ctx, sample_dir):
    """Create a job for a single sample"""
    print("Checking versions of references..")
    fixer = Ref_Updater(ctx.obj['config'], ctx.obj['log'])
    fixer.update_refs()
    print("Version check done. Creating sbatch job")
    worker = Job_Creator(sample_dir, ctx.obj['config'], ctx.obj['log'])
    worker.sample_job()
    done()

@root.group()
@click.pass_context
def scrape(ctx):
  """Parses analysis results and uploads to database"""
  pass

@scrape.command()
@click.argument('sample_dir')
@click.pass_context
def sample(ctx, sample_dir):
  """Parse results from analysing a single sample"""
  garbageman = Scraper(sample_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_sample()
  done()

@scrape.command()
@click.argument('project_dir')
@click.pass_context
def project(ctx, project_dir):
  """Parse results from analysing a single project"""
  garbageman = Scraper(project_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_project()
  done()

@root.command()
@click.pass_context
def update(ctx):
  """Verifies & updates all references"""
  fixer = Ref_Updater(ctx.obj['config'], ctx.obj['log'])
  fixer.update_refs()
  done()

@root.command()
@click.pass_context
def view(ctx):
  """Starts a webserver at http://127.0.0.1:5000/microSALT"""
  app.run()
