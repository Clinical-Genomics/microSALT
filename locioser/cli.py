"""This is the main entry point of locioser. Current commands are analyze and store
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os

from pkg_resources import iter_entry_points
from locioser import __version__

@click.group()
@click.version_option(__version__)
@click.pass_context
def cli(ctx, config_file):
	""" Tool for the Automation of Storage and Analyses """
	ctx.obj = {}
        with open("{}/config.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
          config = yaml.load(conf)

#Add subcommands dynamically to the CLI
for entry_point in iter_entry_points('locioser.subcommands'):
    cli.add_command(entry_point.load())
