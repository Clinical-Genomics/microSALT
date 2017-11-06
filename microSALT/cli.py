"""This is the main entry point of microSALT. Current commands are analyze and store
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pdb
import yaml

from pkg_resources import iter_entry_points
from microSALT import __version__
from microSALT import job_creator
from microSALT import scraper
from microSALT import exporter

@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
    """ Fundamental MLST pipeline """
    with open("{}/config.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
      config = yaml.load(conf)
    ctx.obj = config

@root.command()
@click.argument('indir')
@click.argument('organism')
@click.pass_context
def create_job(ctx, indir, organism):
    boss = job_creator.Job_Creator(indir, organism, ctx.obj)
    boss.create_job()

@root.command()
@click.argument('infile')
@click.pass_context
def scrape(ctx, infile):
  garbageman = scraper.Scraper(infile, ctx.obj)
  garbageman.scrape_blast_loci()

@root.command()
@click.pass_context
def export(ctx):
  ferryman = exporter.Exporter(ctx.obj)
  ferryman.std_export()
