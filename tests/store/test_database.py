import pytest

from unittest.mock import patch
from sqlalchemy import inspect as sa_inspect

from microSALT.exc.exceptions import RefUpdateLockError
from microSALT.store.db_manipulator import DB_Manipulator, _resolve_orm_table
from microSALT.store.orm_models import Reports, Samples, SystemLock


def test_create_every_table(dbm):
    inspector = sa_inspect(dbm.engine)
    assert inspector.has_table("samples")
    assert inspector.has_table("seq_types")
    assert inspector.has_table("resistances")
    assert inspector.has_table("expacs")
    assert inspector.has_table("projects")
    assert inspector.has_table("reports")
    assert inspector.has_table("collections")


def test_add_rec(caplog, profile_dbm):
    dbm = profile_dbm
    # Profile table
    dbm.add_rec(
        {
            "ST": "130",
            "arcC": "6",
            "aroE": "57",
            "glpF": "45",
            "gmk": "2",
            "pta": "7",
            "tpi": "58",
            "yqiL": "52",
        },
        dbm.profiles["staphylococcus_aureus"],
    )
    assert len(dbm.read_records(dbm.profiles["staphylococcus_aureus"], {"ST": "130"})) == 1
    assert len(dbm.read_records(dbm.profiles["staphylococcus_aureus"], {"ST": "-1"})) == 0

    # Novel table
    dbm.add_rec(
        {
            "ST": "130",
            "arcC": "6",
            "aroE": "57",
            "glpF": "45",
            "gmk": "2",
            "pta": "7",
            "tpi": "58",
            "yqiL": "52",
        },
        dbm.novel["staphylococcus_aureus"],
    )
    assert len(dbm.read_records(dbm.novel["staphylococcus_aureus"], {"ST": "130"})) == 1
    assert len(dbm.read_records(dbm.novel["staphylococcus_aureus"], {"ST": "-1"})) == 0

    # ORM tables
    dbm.add_rec({"CG_ID_sample": "ADD1234A1"}, "Samples")
    assert len(dbm.read_records("Samples", {"CG_ID_sample": "ADD1234A1"})) > 0
    assert len(dbm.read_records("Samples", {"CG_ID_sample": "XXX1234A10"})) == 0

    dbm.add_rec({"CG_ID_sample": "ADD1234A1", "loci": "mdh", "contig_name": "NODE_1"}, "Seq_types")
    assert (
        len(
            dbm.read_records(
                "Seq_types", {"CG_ID_sample": "ADD1234A1", "loci": "mdh", "contig_name": "NODE_1"}
            )
        )
        > 0
    )
    assert (
        len(
            dbm.read_records(
                "Seq_types", {"CG_ID_sample": "XXX1234A10", "loci": "mdh", "contig_name": "NODE_1"}
            )
        )
        == 0
    )

    dbm.add_rec(
        {
            "CG_ID_sample": "ADD1234A1",
            "gene": "Type 1",
            "instance": "Type 1",
            "contig_name": "NODE_1",
        },
        "Resistances",
    )
    assert (
        len(
            dbm.read_records(
                "Resistances",
                {
                    "CG_ID_sample": "ADD1234A1",
                    "gene": "Type 1",
                    "instance": "Type 1",
                    "contig_name": "NODE_1",
                },
            )
        )
        > 0
    )
    assert (
        len(
            dbm.read_records(
                "Resistances",
                {
                    "CG_ID_sample": "XXX1234A10",
                    "gene": "Type 1",
                    "instance": "Type 1",
                    "contig_name": "NODE_1",
                },
            )
        )
        == 0
    )

    dbm.add_rec(
        {
            "CG_ID_sample": "ADD1234A1",
            "gene": "Type 1",
            "instance": "Type 1",
            "contig_name": "NODE_1",
        },
        "Expacs",
    )
    assert (
        len(
            dbm.read_records(
                "Expacs",
                {
                    "CG_ID_sample": "ADD1234A1",
                    "gene": "Type 1",
                    "instance": "Type 1",
                    "contig_name": "NODE_1",
                },
            )
        )
        > 0
    )
    assert (
        len(
            dbm.read_records(
                "Expacs",
                {
                    "CG_ID_sample": "XXX1234A10",
                    "gene": "Type 1",
                    "instance": "Type 1",
                    "contig_name": "NODE_1",
                },
            )
        )
        == 0
    )

    dbm.add_rec({"CG_ID_project": "ADD1234"}, "Projects")
    assert len(dbm.read_records("Projects", {"CG_ID_project": "ADD1234"})) > 0
    assert len(dbm.read_records("Projects", {"CG_ID_project": "XXX1234"})) == 0

    dbm.add_rec({"CG_ID_project": "ADD1234", "version": "1"}, "Reports")
    assert len(dbm.read_records("Reports", {"CG_ID_project": "ADD1234", "version": "1"})) > 0
    assert len(dbm.read_records("Reports", {"CG_ID_project": "XXX1234", "version": "1"})) == 0

    dbm.add_rec({"CG_ID_sample": "ADD1234", "ID_collection": "MyCollectionFolder"}, "Collections")
    assert (
        len(
            dbm.read_records(
                "Collections", {"CG_ID_sample": "ADD1234", "ID_collection": "MyCollectionFolder"}
            )
        )
        > 0
    )
    assert (
        len(
            dbm.read_records(
                "Collections", {"CG_ID_sample": "XXX1234", "ID_collection": "MyCollectionFolder"}
            )
        )
        == 0
    )

    caplog.clear()
    dbm.add_rec({"CG_ID_sample": "ADD1234A1"}, "An_entry_that_does_not_exist")
    assert "Attempted to access table" in caplog.text


@patch("sys.exit")
def test_upd_rec(sysexit, caplog, dbm):
    dbm.add_rec({"CG_ID_sample": "UPD1234A1"}, "Samples")
    assert len(dbm.read_records("Samples", {"CG_ID_sample": "UPD1234A1"})) == 1
    assert len(dbm.read_records("Samples", {"CG_ID_sample": "UPD1234A2"})) == 0

    dbm.upd_rec({"CG_ID_sample": "UPD1234A1"}, "Samples", {"CG_ID_sample": "UPD1234A2"})
    assert len(dbm.read_records("Samples", {"CG_ID_sample": "UPD1234A1"})) == 0
    assert len(dbm.read_records("Samples", {"CG_ID_sample": "UPD1234A2"})) == 1

    dbm.upd_rec({"CG_ID_sample": "UPD1234A2"}, "Samples", {"CG_ID_sample": "UPD1234A1"})

    caplog.clear()
    dbm.add_rec({"CG_ID_sample": "UPD1234A1_uniq", "Customer_ID_sample": "cust000"}, "Samples")
    dbm.add_rec({"CG_ID_sample": "UPD1234A2_uniq", "Customer_ID_sample": "cust000"}, "Samples")
    dbm.upd_rec({"Customer_ID_sample": "cust000"}, "Samples", {"Customer_ID_sample": "cust030"})
    dbm.upd_rec({"Customer_ID_sample": "cust000"}, "Samples", {"Customer_ID_sample": "cust030"})
    assert "More than 1 record found" in caplog.text


def test_allele_ranker(profile_dbm, unpack_db_json):
    dbm = profile_dbm
    dbm.add_rec(
        {
            "CG_ID_sample": "MLS1234A1",
            "CG_ID_project": "MLS1234",
            "organism": "staphylococcus_aureus",
        },
        "Samples",
    )
    assert dbm.read_st("MLS1234A1") == 130
    best_alleles = {
        "arcC": {"contig_name": "NODE_1", "allele": 6},
        "aroE": {"contig_name": "NODE_1", "allele": 57},
        "glpF": {"contig_name": "NODE_1", "allele": 45},
        "gmk": {"contig_name": "NODE_1", "allele": 2},
        "pta": {"contig_name": "NODE_1", "allele": 7},
        "tpi": {"contig_name": "NODE_1", "allele": 58},
        "yqiL": {"contig_name": "NODE_1", "allele": 52},
    }
    assert dbm.read_best_alleles("MLS1234A1") == best_alleles

    for entry in unpack_db_json("sampleinfo_mlst.json"):
        entry["allele"] = 0
        entry["CG_ID_sample"] = "MLS1234A2"
        dbm.add_rec(entry, "Seq_types")
    assert dbm.read_st("MLS1234A2") == -1


def test_get_and_set_report(dbm):
    # Clean up any leftover data from prior runs to keep the test idempotent.
    dbm.session.query(Reports).filter(Reports.CG_ID_project == "ADD1234").delete()
    dbm.session.query(Samples).filter(Samples.CG_ID_sample == "ADD1234A1").delete()
    dbm.session.commit()

    dbm.add_rec({"CG_ID_sample": "ADD1234A1", "method_sequencing": "1000:1"}, "Samples")
    dbm.add_rec({"CG_ID_project": "ADD1234", "version": "1"}, "Reports")
    assert dbm.read_report("ADD1234").version == 1

    dbm.upd_rec(
        {"CG_ID_sample": "ADD1234A1", "method_sequencing": "1000:1"},
        "Samples",
        {"CG_ID_sample": "ADD1234A1", "method_sequencing": "1000:2"},
    )
    dbm.set_report("ADD1234")
    assert dbm.read_report("ADD1234").version != 1


@patch("sys.exit")
def test_purge_rec(sysexit, caplog, dbm):
    dbm.add_rec({"CG_ID_sample": "UPD1234A1"}, "Samples")
    dbm.delete_records("UPD1234A1", "Collections")

    caplog.clear()
    dbm.delete_records("UPD1234A1", "Not_Samples_nor_Collections")
    assert "Incorrect type" in caplog.text


def test_top_index(dbm):
    dbm.add_rec({"CG_ID_sample": "Uniq_ID_123", "total_reads": 100}, "Samples")
    dbm.add_rec({"CG_ID_sample": "Uniq_ID_321", "total_reads": 100}, "Samples")
    ti_returned = dbm.read_top_index("Samples", {"total_reads": "100"}, "total_reads")
    assert ti_returned == 100

    ti_missing = dbm.read_top_index("Samples", {"total_reads": "99999"}, "total_reads")
    assert ti_missing == -1


def test_query_rec(dbm):
    dbm.add_rec({"CG_ID_sample": "QRY_001"}, "Samples")
    dbm.add_rec({"CG_ID_sample": "QRY_002"}, "Samples")

    hits = dbm.read_records("Samples", {"CG_ID_sample": "QRY_001"})
    assert len(hits) == 1
    assert hits[0].CG_ID_sample == "QRY_001"

    no_hits = dbm.read_records("Samples", {"CG_ID_sample": "DOES_NOT_EXIST"})
    assert len(no_hits) == 0

    multi_filter = dbm.read_records("Samples", {"CG_ID_sample": "QRY_001", "ST": None})
    assert len(multi_filter) == 1


def test_get_columns(dbm):
    cols = dbm.read_columns("Samples")
    assert isinstance(cols, dict)
    assert "CG_ID_sample" in cols
    assert "organism" in cols


def test_exists(dbm):
    dbm.add_rec({"CG_ID_sample": "EXS_001"}, "Samples")

    assert dbm.exists("Samples", {"CG_ID_sample": "EXS_001"}) is True
    assert dbm.exists("Samples", {"CG_ID_sample": "DOES_NOT_EXIST"}) is False


def test_add_rec_unknown_table(caplog, dbm):
    dbm.add_rec({"CG_ID_sample": "ADD1234A1"}, "An_entry_that_does_not_exist")
    assert "Attempted to access table" in caplog.text


def test_resolve_orm_table_unknown():
    with pytest.raises(KeyError):
        _resolve_orm_table("NonExistentTable")


def test_resolve_orm_table_known():
    assert _resolve_orm_table("Samples") is Samples


def test_populate_profiletable(profile_dbm):
    """populate_profiletable bulk-inserts all data rows from the profile file."""
    # Given: a profile table created from a file with two STs (ST=1 and ST=130).
    dbm = profile_dbm
    table = dbm.profiles["staphylococcus_aureus"]

    # When: populate_profiletable has been called by the profile_dbm fixture.

    # Then: the table contains exactly the two rows from the file.
    rows = dbm.session.query(table).all()
    assert len(rows) == 2

    sts = {row.ST for row in rows}
    assert 130 in sts
    assert 1 in sts

    # Then: allele values for ST=130 match the source file exactly.
    st130 = next(r for r in rows if r.ST == 130)
    assert st130.arcC == 6
    assert st130.aroE == 57
    assert st130.glpF == 45
    assert st130.gmk == 2
    assert st130.pta == 7
    assert st130.tpi == 58
    assert st130.yqiL == 52


def test_refresh_profiletable_same_schema(tmp_profiles_dir, profile_dbm):
    """refresh_profiletable with unchanged columns truncates and reloads data."""
    # Given: a populated profile table (ST=1, ST=130) and a new file with the
    # same column layout but different data (ST=99 only).
    dbm = profile_dbm
    new_content = "ST\tarcC\taroE\tglpF\tgmk\tpta\ttpi\tyqiL\n" "99\t3\t3\t3\t3\t3\t3\t3\n"
    (tmp_profiles_dir / "staphylococcus_aureus").write_text(new_content)

    # When: refresh_profiletable is called with the updated file.
    dbm.refresh_profiletable("staphylococcus_aureus")

    # Then: only the new row is present.
    table = dbm.profiles["staphylococcus_aureus"]
    rows = dbm.session.query(table).all()
    assert len(rows) == 1
    assert rows[0].ST == 99

    # Then: old data has been cleared.
    old_rows = dbm.session.query(table).filter(table.c.ST == 130).all()
    assert len(old_rows) == 0


def test_refresh_profiletable_schema_change(tmp_profiles_dir, profile_dbm):
    """refresh_profiletable drops and rebuilds the table when column layout changes."""
    # Given: a populated profile table with columns ST+7 loci (including yqiL)
    # and a new file that renames yqiL to renamedLocus within the 8-column window.
    dbm = profile_dbm
    new_content = "ST\tarcC\taroE\tglpF\tgmk\tpta\ttpi\trenamedLocus\n" "200\t1\t2\t3\t4\t5\t6\t7\n"
    (tmp_profiles_dir / "staphylococcus_aureus").write_text(new_content)

    # When: refresh_profiletable detects the schema change and does a full drop/recreate.
    dbm.refresh_profiletable("staphylococcus_aureus")

    # Then: the table schema reflects the new column layout.
    table = dbm.profiles["staphylococcus_aureus"]
    cols = list(table.c.keys())
    assert "renamedLocus" in cols
    assert "yqiL" not in cols
    assert len(cols) == 8  # ST + 7 loci

    # Then: the table is populated with data from the new file.
    rows = dbm.session.query(table).all()
    assert len(rows) == 1
    assert rows[0].ST == 200


def test_acquire_ref_lock(clean_lock_dbm):
    """acquire_ref_lock inserts a lock row into the database."""
    dbm = clean_lock_dbm

    # Given: no lock exists.
    assert dbm.session.query(SystemLock).filter_by(lock_name="ref_update").scalar() is None

    # When: the lock is acquired.
    dbm.acquire_ref_lock()

    # Then: a SystemLock row with lock_name='ref_update' is present.
    lock = dbm.session.query(SystemLock).filter_by(lock_name="ref_update").scalar()
    assert lock is not None
    assert lock.lock_name == "ref_update"
    assert lock.acquired_at is not None


def test_acquire_ref_lock_raises_when_already_held(clean_lock_dbm):
    """acquire_ref_lock raises RefUpdateLockError if the lock is already held."""
    dbm = clean_lock_dbm

    # Given: the lock is already held.
    dbm.acquire_ref_lock()

    # When / Then: a second acquire attempt raises immediately.
    with pytest.raises(RefUpdateLockError):
        dbm.acquire_ref_lock()


def test_release_ref_lock(clean_lock_dbm):
    """release_ref_lock removes the lock row so the lock can be re-acquired."""
    dbm = clean_lock_dbm

    # Given: the lock is held.
    dbm.acquire_ref_lock()
    assert dbm.session.query(SystemLock).filter_by(lock_name="ref_update").scalar() is not None

    # When: the lock is released.
    dbm.release_ref_lock()

    # Then: no lock row exists.
    assert dbm.session.query(SystemLock).filter_by(lock_name="ref_update").scalar() is None

    # Then: the lock can be acquired again without raising.
    dbm.acquire_ref_lock()
    assert dbm.session.query(SystemLock).filter_by(lock_name="ref_update").scalar() is not None


def test_check_ref_lock(clean_lock_dbm):
    """check_ref_lock raises when the lock is held and is silent when it is not."""
    dbm = clean_lock_dbm

    # Given: no lock exists.
    # When / Then: check passes silently.
    dbm.check_ref_lock()  # must not raise

    # Given: the lock is acquired.
    dbm.acquire_ref_lock()

    # When / Then: check raises RefUpdateLockError.
    with pytest.raises(RefUpdateLockError):
        dbm.check_ref_lock()

    # Given: the lock is released.
    dbm.release_ref_lock()

    # When / Then: check passes silently again.
    dbm.check_ref_lock()  # must not raise
