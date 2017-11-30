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
