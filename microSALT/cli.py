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
from microSALT.server.views import app
from microSALT.utils.scraper import Scraper
from microSALT.utils.job_creator import Job_Creator
from microSALT.utils.ref_updater import Ref_Updater

if config == '':
  print("ERROR: No properly set-up config under neither ~/.microSALT/config.json nor envvar MICROSALT_CONFIG. Exiting.")
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
  ch.setLevel(logging.INFO)
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
    worker.project_job(single_sample=True)
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
