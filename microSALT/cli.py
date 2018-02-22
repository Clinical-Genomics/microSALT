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
def start(ctx):
  """Starts analysis of project/sample"""
  pass

@start.command()
@click.argument('project_dir')
@click.pass_context
def project(ctx, project_dir):
  """Analyze a whole project"""
  print("Checking versions of references..")
  fixer = Ref_Updater(ctx.obj['config'], ctx.obj['log'])
  fixer.update_refs()
  print("Version check done. Creating sbatch jobs")
  manager = Job_Creator(project_dir, ctx.obj['config'], ctx.obj['log'])
  manager.project_job()
  done() 

@start.command()
@click.argument('sample_dir')
@click.pass_context
def sample(ctx, sample_dir):
    """Analyze a single sample"""
    print("Checking versions of references..")
    fixer = Ref_Updater(ctx.obj['config'], ctx.obj['log'])
    fixer.update_refs()
    print("Version check done. Creating sbatch job")
    worker = Job_Creator(sample_dir, ctx.obj['config'], ctx.obj['log'])
    worker.project_job(single_sample=True)
    done()

@root.group()
@click.pass_context
def finish(ctx):
  """Uploads analysis and generates results"""
  pass

@finish.command()
@click.argument('sample_dir')
@click.pass_context
def sample(ctx, sample_dir):
  """Parse results from analysing a single sample"""
  garbageman = Scraper(sample_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_sample()
  done()

@finish.command()
@click.argument('project_dir')
@click.pass_context
def project(ctx, project_dir):
  """Parse results from analysing a single project"""
  garbageman = Scraper(project_dir, ctx.obj['config'], ctx.obj['log'])
  garbageman.scrape_project()
  done()

@root.command()
@click.argument('project_name')
@click.pass_context
def report(ctx, project_name):
  """Re-generates reports for a project"""
  codemonkey = Reporter(ctx.obj['config'], ctx.obj['log'])
  codemonkey.gen_pdf(project_name)
  codemonkey.gen_csv(project_name)
  done() 
