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
class Seq_types:

  def __init__(self, metadata):

    self.table = Table('seq_types', metadata,
        Column('CG_ID_sample', String(15), ForeignKey(metadata.tables['samples'].c.CG_ID_sample), primary_key=True),
        Column('loci', String(10)),
        Column('assumed_ST', SmallInteger),
        Column('allele', SmallInteger),
        Column('haplotype', String(5)),
        Column('contig_name', String(20), primary_key=True),
        Column('contig_length', Integer),
        Column('contig_coverage', Float(6,2)),
        Column('identity', Float(3,2)),
        Column('evalue', String(10)),
        Column('bitscore', SmallInteger),
        Column('contig_start', Integer),
        Column('contig_end', Integer),
        Column('loci_start', Integer),
        Column('loci_end', Integer),
      )
