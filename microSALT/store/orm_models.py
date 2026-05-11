"""Samples table definition
By: Isak Sylvin, @sylvinite"""

from sqlalchemy import (
    Boolean,
    Column,
    DateTime,
    Float,
    ForeignKey,
    Integer,
    SmallInteger,
    String,
)
from sqlalchemy.orm import declarative_base, relationship

Base = declarative_base()


class Samples(Base):
    __tablename__ = "samples"
    seq_types = relationship("Seq_types", back_populates="samples")
    projects = relationship("Projects", back_populates="samples")
    resistances = relationship("Resistances", back_populates="samples")
    # steps = relationship("Steps", back_populates="samples")
    expacs = relationship("Expacs", back_populates="samples")

    CG_ID_sample = Column(String(32), primary_key=True, nullable=False)
    CG_ID_project = Column(String(32), ForeignKey("projects.CG_ID_project"))
    Customer_ID_sample = Column(String(128))
    organism = Column(String(100))
    ST = Column(SmallInteger, default=-1)
    pubmlst_ST = Column(SmallInteger, default=-1)
    date_analysis = Column(DateTime)
    genome_length = Column(Integer, default=-1)
    gc_percentage = Column(Float(3, 2), default=0.0)
    n50 = Column(Integer, default=-1)
    contigs = Column(Integer, default=-1)
    priority = Column(String(20))

    total_reads = Column(Integer)  # Fetch from bcl2fastq
    insert_size = Column(Integer)
    duplication_rate = Column(Float)
    mapped_rate = Column(Float)
    coverage_10x = Column(Float)
    coverage_30x = Column(Float)
    coverage_50x = Column(Float)
    coverage_100x = Column(Float)
    average_coverage = Column(Float)
    reference_genome = Column(String(32))
    reference_length = Column(Integer)

    application_tag = Column(String(15))
    date_arrival = Column(DateTime)
    date_analysis = Column(DateTime)
    date_sequencing = Column(DateTime)
    date_libprep = Column(DateTime)
    method_sequencing = Column(String(128))
    method_libprep = Column(String(128))


class Seq_types(Base):
    __tablename__ = "seq_types"
    samples = relationship("Samples", back_populates="seq_types")

    CG_ID_sample = Column(String(32), ForeignKey("samples.CG_ID_sample"), primary_key=True)
    loci = Column(String(10), primary_key=True)
    allele = Column(SmallInteger)
    contig_name = Column(String(20), primary_key=True)
    contig_length = Column(Integer)
    contig_coverage = Column(Float(6, 2))
    identity = Column(Float(3, 2), default=0.0)
    span = Column(Float(3, 2), default=0.0)
    evalue = Column(String(10))
    bitscore = Column(SmallInteger)
    subject_length = Column(Integer)
    st_predictor = Column(Boolean, default=0)
    contig_start = Column(Integer)
    contig_end = Column(Integer)


class Resistances(Base):
    __tablename__ = "resistances"
    samples = relationship("Samples", back_populates="resistances")

    CG_ID_sample = Column(String(32), ForeignKey("samples.CG_ID_sample"), primary_key=True)
    gene = Column(String(50), primary_key=True)
    instance = Column(String(120), primary_key=True)
    contig_name = Column(String(20), primary_key=True)
    contig_length = Column(Integer)
    contig_coverage = Column(Float(6, 2))
    identity = Column(Float(3, 2), default=0.0)
    span = Column(Float(3, 2), default=0.0)
    evalue = Column(String(10))
    bitscore = Column(SmallInteger)
    subject_length = Column(Integer)
    reference = Column(String(40))
    resistance = Column(String(120))
    contig_start = Column(Integer)
    contig_end = Column(Integer)


class Expacs(Base):
    __tablename__ = "expacs"
    samples = relationship("Samples", back_populates="expacs")

    CG_ID_sample = Column(String(32), ForeignKey("samples.CG_ID_sample"), primary_key=True)
    gene = Column(String(50), primary_key=True)
    instance = Column(String(120), primary_key=True)
    contig_name = Column(String(20), primary_key=True)
    contig_length = Column(Integer)
    contig_coverage = Column(Float(6, 2))
    identity = Column(Float(3, 2), default=0.0)
    span = Column(Float(3, 2), default=0.0)
    evalue = Column(String(10))
    bitscore = Column(SmallInteger)
    subject_length = Column(Integer)
    reference = Column(String(40))
    virulence = Column(String(120))
    contig_start = Column(Integer)
    contig_end = Column(Integer)


class Projects(Base):
    __tablename__ = "projects"
    samples = relationship("Samples", back_populates="projects")
    reports = relationship("Reports", back_populates="projects")

    CG_ID_project = Column(String(32), primary_key=True, nullable=False)
    Customer_ID_project = Column(String(32))
    Customer_ID = Column(String(32))


class Versions(Base):
    __tablename__ = "versions"

    name = Column(String(45), primary_key=True, nullable=False)
    version = Column(String(10))


# Keeps and aggregate step string, makes a new version whenever one is not found
class Reports(Base):
    __tablename__ = "reports"
    projects = relationship("Projects", back_populates="reports")

    CG_ID_project = Column(String(32), ForeignKey("projects.CG_ID_project"), primary_key=True)
    steps_aggregate = Column(String(100))
    date = Column(DateTime)
    version = Column(Integer, default=1, primary_key=True)


class Collections(Base):
    __tablename__ = "collections"

    ID_collection = Column(String(32), primary_key=True)
    CG_ID_sample = Column(String(32), primary_key=True)


class SystemLock(Base):
    """Rows in this table act as advisory locks for long-running operations.

    A row with lock_name='ref_update' signals that a reference update is in
    progress.  No other processes should modify or read profile tables while
    this lock is held.
    """

    __tablename__ = "system_locks"

    lock_name = Column(String(60), primary_key=True, nullable=False)
    acquired_at = Column(DateTime, nullable=False)
