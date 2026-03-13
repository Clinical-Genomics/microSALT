import copy

import pytest
from sqlalchemy import inspect as sa_inspect

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.store.orm_models import SystemLock


@pytest.fixture
def tmp_profiles_dir(tmp_path):
    """Creates a temporary profiles directory with a staphylococcus_aureus profile file.

    The file contains two rows: ST=130 (matching sampleinfo_mlst.json alleles) and ST=1.
    """
    content = (
        "ST\tarcC\taroE\tglpF\tgmk\tpta\ttpi\tyqiL\n"
        "1\t1\t1\t1\t1\t1\t1\t1\n"
        "130\t6\t57\t45\t2\t7\t58\t52\n"
    )
    (tmp_path / "staphylococcus_aureus").write_text(content)
    return tmp_path


@pytest.fixture
def profile_dbm(config, logger, tmp_profiles_dir, unpack_db_json):
    """DB_Manipulator with profile/novel tables freshly built from tmp_profiles_dir.

    Uses the shared SQLite database but drops and recreates all profile/novel tables
    on each invocation so tests start with known, clean data.
    """
    cfg = copy.deepcopy(config)
    cfg.folders.profiles = str(tmp_profiles_dir)

    dbm = DB_Manipulator(log=logger, folders=cfg.folders, threshold=cfg.threshold)
    dbm.create_tables()

    inspector = sa_inspect(dbm.engine)

    # Drop and recreate profile tables with data from tmp_profiles_dir.
    for name, table in list(dbm.profiles.items()):
        if inspector.has_table(table.name):
            table.drop(dbm.engine)
        table.create(dbm.engine)
        dbm.populate_profiletable(name, table)

    # Drop and recreate novel tables (empty — entries are written by application logic).
    for name, table in list(dbm.novel.items()):
        if inspector.has_table(table.name):
            table.drop(dbm.engine)
        table.create(dbm.engine)

    for entry in unpack_db_json("sampleinfo_projects.json"):
        dbm.add_to_session(dbm.add_project(**entry))
    dbm.commit_session()
    for entry in unpack_db_json("sampleinfo_mlst.json"):
        dbm.add_to_session(dbm.add_seq_type(**entry))
    dbm.commit_session()

    return dbm


@pytest.fixture
def clean_lock_dbm(dbm):
    """Wraps dbm and guarantees the ref_update lock is absent before and after each test.

    Without this, a test that fails mid-way while holding the lock would leave dirty
    state in the shared SQLite DB and cause unrelated tests to fail.
    """
    dbm.session.query(SystemLock).filter_by(lock_name="ref_update").delete()
    dbm.session.commit()
    yield dbm
    dbm.session.query(SystemLock).filter_by(lock_name="ref_update").delete()
    dbm.session.commit()
