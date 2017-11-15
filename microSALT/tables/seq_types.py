"""Sequencing types (blast results) table definitions
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

from sqlalchemy import *
from microSALT import Base
from sqlalchemy.orm import relationship

class Seq_types(Base):
  __tablename__ = 'seq_types'
  CG_ID_sample = Column(String(15), ForeignKey('samples.CG_ID_sample'), primary_key=True)
  samples = relationship('Samples', back_populates='seq_types')
  loci = Column(String(10))
  allele = Column(SmallInteger)
  haplotype = Column(String(5))
  contig_name = Column(String(20), primary_key=True)
  contig_length = Column(Integer)
  contig_coverage = Column(Float(6,2))
  identity = Column(Float(3,2))
  evalue = Column(String(10))
  bitscore = Column(SmallInteger)
  contig_start = Column(Integer)
  contig_end = Column(Integer)
  loci_start = Column(Integer)
  loci_end = Column(Integer)
