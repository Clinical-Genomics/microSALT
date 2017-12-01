"""Samples table definition
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

from sqlalchemy import *
from sqlalchemy.orm import relationship
from microSALT import db

class Samples(db.Model):

   __tablename__ = 'samples'
   CG_ID_sample = db.Column(db.String(15), primary_key=True, nullable=False)
   seq_types = relationship("Seq_types", back_populates="samples")
   CG_ID_project = db.Column(db.String(15))
   Customer_ID_sample = db.Column(db.String(15))
   Customer_ID_project = db.Column(db.String(15))
   date_ordered = db.Column(db.DateTime)
   date_qc = db.Column(db.DateTime)
   date_analysis = db.Column(db.DateTime)
   organism = db.Column(db.String(30))
   ST = db.Column(db.SmallInteger, default=-1)

class Seq_types(db.Model):
  __tablename__ = 'seq_types'
  CG_ID_sample = db.Column(db.String(15), ForeignKey('samples.CG_ID_sample'), primary_key=True)
  samples = relationship('Samples', back_populates='seq_types')
  loci = db.Column(db.String(10), primary_key=True)
  allele = db.Column(db.SmallInteger)
  haplotype = db.Column(db.String(5))
  contig_name = db.Column(db.String(20), primary_key=True)
  contig_length = db.Column(db.Integer)
  contig_coverage = db.Column(db.Float(6,2))
  identity = db.Column(db.Float(3,2))
  evalue = db.Column(db.String(10))
  bitscore = db.Column(db.SmallInteger)
  contig_start = db.Column(db.Integer)
  contig_end = db.Column(db.Integer)
  loci_start = db.Column(db.Integer)
  loci_end = db.Column(db.Integer)

