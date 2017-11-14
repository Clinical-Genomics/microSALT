"""Samples table definition
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

from sqlalchemy import *
from microSALT import Base
from sqlalchemy.orm import relationship


class Samples(Base):

   __tablename__ = 'samples'
   CG_ID_sample = Column(String(15), primary_key=True, nullable=False)
   seq_types = relationship("Seq_types", back_populates="samples")
   CG_ID_project = Column(String(15))
   Customer_ID_sample = Column(String(15))
   Customer_ID_project = Column(String(15))
   date_analysis = Column(DateTime)
   organism = Column(String(30))
   ST = Column(SmallInteger, default=-1)
