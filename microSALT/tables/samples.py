""" Initial script to deliver and fetch data from a database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import re
from sqlalchemy import *
import yaml

import pdb # debug

#TODO: Rewrite all pushes/queries through session+commit
class Samples:
  def __init__(self, metadata):

    self.table = Table('samples', metadata,
        Column('CG_ID_sample', String(15), primary_key=True, nullable=False),
        Column('CG_ID_project', String(15)),
        Column('Customer_ID_sample', String(15)),
        Column('Customer_ID_project', String(15)),
        Column('date_analysis', DateTime),
        Column('organism', String(30)),
        Column('ST', SmallInteger, default=-1),
      )
