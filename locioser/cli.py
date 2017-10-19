"""This is the main entry point of locioser. Current commands are analyze and store
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pdb
import yaml

from pkg_resources import iter_entry_points
from locioser import __version__
from locioser import job_creator

@click.group()
@click.version_option(__version__)
@click.pass_context
def root(ctx):
    """ Fundamental MLST pipeline """
    ctx.obj = {}

@root.command()
@click.argument('indir')
def create_job(indir):
    boss = job_creator.Job_Creator(indir)
    boss.create_job()
