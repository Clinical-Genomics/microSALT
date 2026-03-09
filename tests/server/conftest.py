import pytest
from datetime import datetime

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

import microSALT.store.database as db_module
from microSALT.store.orm_models import (
    Base,
    Projects,
    Reports,
    Samples,
    Seq_types,
    Versions,
)


@pytest.fixture(scope="module")
def db_engine():
    """In-memory SQLite engine shared for the whole module."""
    engine = create_engine("sqlite:///:memory:")
    Base.metadata.create_all(engine)
    return engine


@pytest.fixture(scope="module")
def db_session(config, db_engine):
    """Scoped session bound to the in-memory engine.

    Depends on ``config`` to ensure initialize_database() (called by that
    fixture) runs *before* we wire db_module.session to the in-memory engine,
    so the in-memory session is not overwritten afterwards.

    Also wires up the module-level session used by get_session() in views.py.
    """
    session_factory = sessionmaker(db_engine)
    session = scoped_session(session_factory)

    # Wire up the database module so get_session() returns this session
    db_module.engine = db_engine
    db_module.session = session

    yield session

    session.remove()


@pytest.fixture(scope="module")
def populated_db(db_session):
    """Populate the in-memory database with a minimal set of test data."""
    # Project
    project = Projects(
        CG_ID_project="AAA1234",
        Customer_ID_project="999999",
        date_ordered=datetime(2020, 7, 14, 15, 3, 51),
        Customer_ID="cust000",
    )
    db_session.add(project)

    # Samples
    sample1 = Samples(
        CG_ID_sample="AAA1234A1",
        CG_ID_project="AAA1234",
        Customer_ID_sample="XXX0000Y1",
        organism="staphylococcus_aureus",
        ST=8,
        pubmlst_ST=-1,
        total_reads=500000,
        insert_size=200,
        duplication_rate=0.05,
        mapped_rate=0.95,
        average_coverage=120.0,
        coverage_10x=0.99,
        coverage_30x=0.97,
        coverage_50x=0.90,
        coverage_100x=0.50,
        genome_length=2800000,
        gc_percentage=32.5,
        n50=250000,
        contigs=12,
    )
    sample2 = Samples(
        CG_ID_sample="AAA1234A2",
        CG_ID_project="AAA1234",
        Customer_ID_sample="XXX0000Y2",
        organism="escherichia_coli",
        ST=131,
        pubmlst_ST=-1,
        total_reads=600000,
        insert_size=180,
        duplication_rate=0.04,
        mapped_rate=0.97,
        average_coverage=150.0,
        coverage_10x=0.99,
        coverage_30x=0.98,
        coverage_50x=0.95,
        coverage_100x=0.70,
        genome_length=5000000,
        gc_percentage=50.0,
        n50=400000,
        contigs=10,
    )
    db_session.add(sample1)
    db_session.add(sample2)

    # Seq types for sample1
    for loci, allele in [("arcC", 6), ("aroE", 57), ("glpF", 45)]:
        st = Seq_types(
            CG_ID_sample="AAA1234A1",
            loci=loci,
            allele=allele,
            contig_name="NODE_1",
            contig_length=592262,
            contig_coverage=150.0,
            identity=100.0,
            span=1.0,
            evalue="0.0",
            bitscore=900,
            subject_length=456,
            st_predictor=True,
            contig_start=100,
            contig_end=556,
        )
        db_session.add(st)

    # Report
    report = Reports(
        CG_ID_project="AAA1234",
        steps_aggregate="abc123",
        date=datetime(2020, 7, 14, 15, 3, 51),
        version=1,
    )
    db_session.add(report)

    # Versions
    v = Versions(name="software_microSALT", version="4.3.0")
    db_session.add(v)

    db_session.commit()

    yield db_session
