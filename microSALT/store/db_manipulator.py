"""Delivers and fetches data from the database
By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import hashlib
import sys
import warnings
from collections import OrderedDict
from datetime import datetime, timezone

from dateutil.parser import parse
from sqlalchemy import DateTime as SADateTime
from sqlalchemy import MetaData, Table, and_, desc, or_, text
from sqlalchemy import inspect as sa_inspect

from microSALT import __version__
from microSALT.config import Folders, Threshold
from microSALT.exc.exceptions import RefUpdateLockError
from microSALT.store.database import get_engine, get_session
from microSALT.store.models import ProfileTable
from microSALT.store.orm_models import (
    Collections,
    Expacs,
    Projects,
    Reports,
    Resistances,
    Samples,
    Seq_types,
    SystemLock,
    Versions,
)

# Maps string table names (as passed by callers) to ORM classes.
_ORM_TABLES = {
    "Collections": Collections,
    "Expacs": Expacs,
    "Projects": Projects,
    "Reports": Reports,
    "Resistances": Resistances,
    "Samples": Samples,
    "Seq_types": Seq_types,
    "Versions": Versions,
}


def _resolve_orm_table(tablename: str):
    """Return the ORM class for *tablename*, raising KeyError on unknown names."""
    if tablename not in _ORM_TABLES:
        raise KeyError(f"Unknown ORM table: {tablename!r}")
    return _ORM_TABLES[tablename]


class DB_Manipulator:
    def __init__(self, log, folders: Folders, threshold: Threshold):
        self.folders = folders
        self.threshold = threshold
        self.logger = log
        self.session = get_session()
        self.engine = get_engine()
        self.metadata = MetaData()
        self.profiles = ProfileTable(
            "profile_", self.metadata, self.folders.profiles, self.logger
        ).tables
        self.novel = ProfileTable(
            "novel_", self.metadata, self.folders.profiles, self.logger
        ).tables
        # Turns off pymysql deprecation warnings until they can update their code
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.create_tables()

    def create_tables(self):
        """Creates all tables individually. A bit more control than usual"""
        inspector = sa_inspect(self.engine)
        if not inspector.has_table("projects"):
            Projects.__table__.create(self.engine)
            self.logger.info("Created projects table")
        if not inspector.has_table("samples"):
            Samples.__table__.create(self.engine)
            self.logger.info("Created samples table")
        if not inspector.has_table("versions"):
            Versions.__table__.create(self.engine)
            self.logger.info("Created versions table")
        if not inspector.has_table("seq_types"):
            Seq_types.__table__.create(self.engine)
            self.logger.info("Created sequencing types table")
        if not inspector.has_table("resistances"):
            Resistances.__table__.create(self.engine)
            self.logger.info("Created resistance table")
        if not inspector.has_table("reports"):
            Reports.__table__.create(self.engine)
            self.logger.info("Created reports table")
        if not inspector.has_table("collections"):
            Collections.__table__.create(self.engine)
            self.logger.info("Created collections table")
        if not inspector.has_table("expacs"):
            Expacs.__table__.create(self.engine)
            self.logger.info("Created ExPEC table")
        if not inspector.has_table("system_locks"):
            SystemLock.__table__.create(self.engine)
            self.logger.info("Created system_locks table")
        for k, v in self.profiles.items():
            if not inspector.has_table(f"profile_{k}"):
                self.profiles[k].create(self.engine)
                self.populate_profiletable(k, v)
                self.add_to_session(self.add_version(name=f"profile_{k}", version="0"))
                self.commit_session()
                self.logger.info(f"Profile table profile_{k} created and populated")
        for k, v in self.novel.items():
            if not inspector.has_table(f"novel_{k}"):
                self.novel[k].create(self.engine)
                self.add_to_session(self.add_version(name=f"novel_{k}", version="0"))
                self.commit_session()
                self.logger.info(f"Profile table novel_{k} initialized")

    def acquire_ref_lock(self):
        """Acquire the reference-update exclusive lock.

        Raises RefUpdateLockError if the lock is already held by another process.
        """
        existing = self.session.query(SystemLock).filter_by(lock_name="ref_update").scalar()
        if existing:
            raise RefUpdateLockError(
                "A reference update is already in progress (lock acquired at {}). "
                "Please try again later.".format(existing.acquired_at)
            )
        self.session.add(SystemLock(lock_name="ref_update", acquired_at=datetime.now(timezone.utc)))
        self.session.commit()
        self.logger.info("Reference update lock acquired")

    def release_ref_lock(self):
        """Release the reference-update exclusive lock."""
        self.session.query(SystemLock).filter_by(lock_name="ref_update").delete()
        self.session.commit()
        self.logger.info("Reference update lock released")

    def check_ref_lock(self):
        """Raise RefUpdateLockError if a reference update is currently in progress."""
        lock = self.session.query(SystemLock).filter_by(lock_name="ref_update").scalar()
        if lock:
            raise RefUpdateLockError(
                "The reference database is currently being updated (started at {}). "
                "Please try again later.".format(lock.acquired_at)
            )

    def add_rec(self, data_dict: dict, tablename) -> None:
        """Adds a record to a non-ORM (ProfileTable) table via a raw Table object."""
        table = tablename
        pk_list = table.primary_key.columns.keys()
        filter_clauses = [table.c[pk] == data_dict[pk] for pk in pk_list]
        exist = self.session.query(table).filter(or_(*filter_clauses)).all()
        if len(exist) == 0:
            data = table.insert()
            for k, v in data_dict.items():
                if isinstance(v, str):
                    try:
                        parse(v, fuzzy=False)
                        data_dict[k] = datetime.strptime(v, "%Y-%m-%d %H:%M:%S")
                    except ValueError as ve:
                        if len(ve.args) > 0 and ve.args[0].startswith("unconverted data remains: "):
                            data_dict[k] = datetime.strptime(v, "%Y-%m-%d %H:%M:%S.%f")
                        else:
                            pass
            self.session.execute(data, data_dict)
            self.session.commit()
            self.logger.info(f"Added entry to table {tablename.fullname}")

    # ------------------------------------------------------------------
    # Per-model factory methods — construct and return an ORM object.
    # Callers are responsible for add_to_session() and commit_session().
    # ------------------------------------------------------------------

    def add_sample(self, **kwargs) -> Samples:
        return Samples(**kwargs)

    def add_project(self, **kwargs) -> Projects:
        return Projects(**kwargs)

    def add_seq_type(self, **kwargs) -> Seq_types:
        return Seq_types(**kwargs)

    def add_resistance(self, **kwargs) -> Resistances:
        return Resistances(**kwargs)

    def add_expac(self, **kwargs) -> Expacs:
        return Expacs(**kwargs)

    def add_report(self, **kwargs) -> Reports:
        return Reports(**kwargs)

    def add_collection(self, **kwargs) -> Collections:
        return Collections(**kwargs)

    def add_version(self, **kwargs) -> Versions:
        return Versions(**kwargs)

    def add_to_session(self, obj) -> None:
        """Coerce string DateTime fields, then stage the object.

        If an object with the same primary key already exists in the database
        the insert is silently skipped (same semantics as the previous
        _add_orm_record behaviour).
        """
        for col in obj.__table__.columns:
            if isinstance(col.type, SADateTime):
                val = getattr(obj, col.name)
                if isinstance(val, str):
                    try:
                        setattr(obj, col.name, datetime.strptime(val, "%Y-%m-%d %H:%M:%S"))
                    except ValueError:
                        setattr(obj, col.name, datetime.strptime(val, "%Y-%m-%d %H:%M:%S.%f"))
        pk_cols = list(obj.__table__.primary_key.columns.keys())
        pk_vals = [getattr(obj, c) for c in pk_cols]
        if None not in pk_vals:
            existing = self.session.get(type(obj), pk_vals if len(pk_vals) > 1 else pk_vals[0])
            if existing is not None:
                return
        self.session.add(obj)

    def commit_session(self) -> None:
        self.session.commit()

    # ------------------------------------------------------------------
    # Per-model update methods
    # ------------------------------------------------------------------

    def update_sample(self, req_dict: dict, upd_dict: dict) -> None:
        """Update a Samples row. Cascades CG_ID_sample renames to child tables."""
        filter_clauses = [getattr(Samples, k) == v for k, v in req_dict.items() if v is not None]
        query = self.session.query(Samples).filter(and_(*filter_clauses))
        if len(query.all()) > 1:
            self.logger.error("More than 1 Samples record found when updating. Exited.")
            sys.exit()
        if "CG_ID_sample" in upd_dict:
            old_id = req_dict.get("CG_ID_sample")
            new_id = upd_dict["CG_ID_sample"]
            if old_id and old_id != new_id:
                for child_table in (Seq_types, Resistances, Expacs, Collections):
                    self.session.query(child_table).filter(
                        child_table.CG_ID_sample == old_id
                    ).update({"CG_ID_sample": new_id})
        query.update(upd_dict)
        self.session.commit()
        self.logger.debug(f"Updated Samples for {req_dict} with {upd_dict}")

    def update_project(self, req_dict: dict, upd_dict: dict) -> None:
        """Update a Projects row."""
        filter_clauses = [getattr(Projects, k) == v for k, v in req_dict.items() if v is not None]
        query = self.session.query(Projects).filter(and_(*filter_clauses))
        if len(query.all()) > 1:
            self.logger.error("More than 1 Projects record found when updating. Exited.")
            sys.exit()
        query.update(upd_dict)
        self.session.commit()
        self.logger.debug(f"Updated Projects for {req_dict} with {upd_dict}")

    def update_version(self, req_dict: dict, upd_dict: dict) -> None:
        """Update a Versions row."""
        filter_clauses = [getattr(Versions, k) == v for k, v in req_dict.items() if v is not None]
        self.session.query(Versions).filter(and_(*filter_clauses)).update(upd_dict)
        self.session.commit()
        self.logger.debug(f"Updated Versions for {req_dict} with {upd_dict}")

    # ------------------------------------------------------------------
    # Per-model delete methods
    # ------------------------------------------------------------------

    def delete_sample(self, cg_id: str) -> None:
        """Delete a sample and all its child rows (seq_types, resistances, expacs)."""
        for obj in self.session.query(Expacs).filter(Expacs.CG_ID_sample == cg_id).all():
            self.session.delete(obj)
        for obj in self.session.query(Seq_types).filter(Seq_types.CG_ID_sample == cg_id).all():
            self.session.delete(obj)
        for obj in self.session.query(Resistances).filter(Resistances.CG_ID_sample == cg_id).all():
            self.session.delete(obj)
        for obj in self.session.query(Samples).filter(Samples.CG_ID_sample == cg_id).all():
            self.session.delete(obj)
        self.session.commit()
        self.logger.info(f"Removed sample {cg_id} and its child rows")

    def delete_sample_results(self, cg_id: str) -> None:
        """Delete only the analysis result rows for a sample (seq_types, resistances, expacs)
        without removing the Samples row itself."""
        for obj in self.session.query(Expacs).filter(Expacs.CG_ID_sample == cg_id).all():
            self.session.delete(obj)
        for obj in self.session.query(Seq_types).filter(Seq_types.CG_ID_sample == cg_id).all():
            self.session.delete(obj)
        for obj in self.session.query(Resistances).filter(Resistances.CG_ID_sample == cg_id).all():
            self.session.delete(obj)
        self.session.commit()
        self.logger.info(f"Cleared analysis results for sample {cg_id}")

    def delete_project(self, name: str) -> None:
        """Delete all samples (and their child rows) belonging to a project."""
        for obj in self.session.query(Expacs).filter(Expacs.CG_ID_sample.like(f"{name}%")).all():
            self.session.delete(obj)
        for obj in (
            self.session.query(Seq_types).filter(Seq_types.CG_ID_sample.like(f"{name}%")).all()
        ):
            self.session.delete(obj)
        for obj in (
            self.session.query(Resistances).filter(Resistances.CG_ID_sample.like(f"{name}%")).all()
        ):
            self.session.delete(obj)
        for obj in self.session.query(Samples).filter(Samples.CG_ID_sample.like(f"{name}%")).all():
            self.session.delete(obj)
        self.session.commit()
        self.logger.info(f"Removed all samples for project {name}")

    def delete_collection(self, name: str) -> None:
        """Delete all entries for the given collection ID."""
        for obj in self.session.query(Collections).filter(Collections.ID_collection == name).all():
            self.session.delete(obj)
        self.session.commit()
        self.logger.info(f"Removed collection {name}")

    def read_records(self, tablename: str, filters: dict[str, str]):
        """Fetches records table, using a primary-key dict with columns as keys.
        Non-PK are ignored"""
        # Non-orm
        if not isinstance(tablename, str):
            table = tablename
            filter_clauses = [table.c[k] == v for k, v in filters.items()]
            return self.session.query(table).filter(or_(*filter_clauses)).all()
        # ORM
        else:
            table = _resolve_orm_table(tablename)
            filter_clauses = [getattr(table, k) == v for k, v in filters.items() if v is not None]
            return self.session.query(table).filter(and_(*filter_clauses)).all()

    def read_top_index(self, table_str: str, filters: dict[str, str], column: str):
        """Fetches the top index from column of table, by applying a dict with columns as keys."""
        table = _resolve_orm_table(table_str)
        filter_clauses = [getattr(table, k) == v for k, v in filters.items() if v is not None]
        entry = (
            self.session.query(table)
            .filter(and_(*filter_clauses))
            .order_by(desc(getattr(table, column)))
            .limit(1)
            .all()
        )
        if entry == []:
            return int(-1)
        else:
            return getattr(entry[0], column)

    def reload_profiletable(self, organism: str):
        """Drop the named profile table, rebuild schema from disk, and reload with fresh data.

        The Python Table object is rebuilt from the current file on disk before the DB
        table is recreated, so schema changes (new or renamed loci columns) are picked up.
        """
        self.logger.debug(f"Reloading profile table for {organism}")
        self.profiles[organism].drop(self.engine)
        self.logger.debug(f"Dropped profile table for {organism}")
        # Rebuild the Table object from the file currently on disk (schema may have changed).
        fresh_metadata = MetaData()
        fresh = ProfileTable("profile_", fresh_metadata, self.folders.profiles, self.logger).tables
        if organism in fresh:
            self.profiles[organism] = fresh[organism]
        self.profiles[organism].create(self.engine)
        self.logger.debug(f"Recreated profile table for {organism}")
        self.populate_profiletable(organism, self.profiles[organism])
        self.logger.debug(f"Populated profile table for {organism}")

    def refresh_profiletable(self, organism: str):
        """Reload profile table content without dropping the table when possible.

        Reads the downloaded CSV header and compares it against the current
        table's columns (first 8). If the schema is unchanged, the table is
        truncated and reloaded in place. If the loci scheme has changed (new
        or renamed columns) the method falls back to a full drop/recreate via
        reload_profiletable() so the schema stays in sync with the CSV.
        """
        table = self.profiles[organism]
        file_path = f"{self.folders.profiles}/{organism}"

        with open(file_path, "r") as fh:
            csv_cols = fh.readline().rstrip().split("\t")[:8]

        current_cols = list(table.c.keys())

        if csv_cols == current_cols:
            self.logger.info(
                f"Schema unchanged for {organism}, truncating and reloading profile table"
            )
            self.session.execute(table.delete())
            self.session.commit()
            self.populate_profiletable(organism, table)
        else:
            self.logger.info(
                f"Schema changed for {organism} ({current_cols} -> {csv_cols}), "
                f"dropping and recreating profile table"
            )
            self.reload_profiletable(organism)

    def populate_profiletable(self, filename: str, table) -> None:
        """Bulk-inserts all data rows from a profile file into an already-created *table*."""
        file_path = f"{self.folders.profiles}/{filename}"
        self.logger.debug(f"Opening profile file: {file_path}")
        keys = list(table.c.keys())
        rows = []
        with open(file_path, "r") as fh:
            head = fh.readline().rstrip().split("\t")
            for raw_line in fh:
                values = raw_line.rstrip().split("\t")
                row = {col: None for col in keys}
                for i, val in enumerate(values[: len(keys)]):
                    row[head[i]] = val
                rows.append(row)
        if not rows:
            self.logger.warning(f"No data rows found in profile file {filename}")
            return
        try:
            self.session.execute(table.insert(), rows)
            self.session.commit()
            self.logger.debug(f"Inserted {len(rows)} rows into profile table for {filename}")
        except Exception as e:
            self.session.rollback()
            self.logger.error(f"Failed to bulk-insert profile data for {filename}: {e}")

    def read_columns(self, tablename: str):
        """Returns all records for a given ORM table"""
        table = _resolve_orm_table(tablename)
        return dict.fromkeys(table.__table__.columns.keys())

    def get_projects_by_cg_id_project(self, cg_id_project_name: str) -> Projects | None:
        """Fetch a Projects record by CG_ID_project."""
        return (
            self.session.query(Projects)
            .filter(Projects.CG_ID_project == cg_id_project_name)
            .scalar()
        )

    def get_collection_by_id(self, collection_id: str) -> Collections | None:
        return (
            self.session.query(Collections)
            .filter(Collections.ID_collection == collection_id)
            .scalar()
        )

    def get_sample_by_cg_id_sample(self, cg_id_sample: str) -> Samples | None:
        return self.session.query(Samples).filter(Samples.CG_ID_sample == cg_id_sample).scalar()

    def read_version(self, name: str):
        """Gets the version from a given name. Should be generalized to return any value for any input"""
        version: Versions | None = (
            self.session.query(Versions).filter(Versions.name == name).scalar()
        )
        if version is None:
            return "0"
        else:
            return version.version

    def read_report(self, name: str) -> Reports | None:
        # Sort based on version
        prev_report: Reports | None = None
        prev_reports = (
            self.session.query(Reports)
            .filter(Reports.CG_ID_project == name)
            .order_by(desc(Reports.version))
            .all()
        )
        if len(prev_reports) > 0:
            prev_report = prev_reports[0]
        return prev_report

    def set_report(self, name: str):
        # Generate string
        totalstring: list[str] = []
        dt = datetime.now()
        default_method = "Not in LIMS"
        samples = (
            self.session.query(Samples)
            .filter(Samples.CG_ID_project == name)
            .order_by(desc(Samples.CG_ID_sample))
            .all()
        )
        for sample in samples:
            if sample.date_libprep:
                totalstring.append(
                    str(datetime.timestamp(sample.date_libprep.replace(tzinfo=timezone.utc)))
                )
            else:
                totalstring.append(
                    str(datetime.timestamp(datetime.min.replace(tzinfo=timezone.utc)))
                )

            if sample.method_libprep:
                totalstring.append(sample.method_libprep)
            else:
                totalstring.append(default_method)

            if sample.date_sequencing:
                totalstring.append(
                    str(datetime.timestamp(sample.date_sequencing.replace(tzinfo=timezone.utc)))
                )
            else:
                totalstring.append(
                    str(datetime.timestamp(datetime.min.replace(tzinfo=timezone.utc)))
                )

            if sample.method_sequencing:
                totalstring.append(sample.method_sequencing)
            else:
                totalstring.append(default_method)

        totalstring.append(__version__)
        totalstring = "".join(totalstring).encode()
        hashstring = hashlib.md5(totalstring).hexdigest()

        if prev_report := self.read_report(name):
            if "steps_aggregate" in dir(prev_report) and prev_report.steps_aggregate != hashstring:
                self.add_to_session(
                    self.add_report(
                        CG_ID_project=name,
                        steps_aggregate=hashstring,
                        date=dt,
                        version=prev_report.version + 1,
                    )
                )
                self.commit_session()
        else:
            self.add_to_session(
                self.add_report(
                    CG_ID_project=name,
                    steps_aggregate=hashstring,
                    date=dt,
                    version=1,
                )
            )
            self.commit_session()

    def set_novel_st(self, overwrite=False, sample=""):
        """Looks at each novel table. See if any record has a profile match in the profile table.
        Updates these based on parameters"""
        prequery = self.session.query(Samples)

        for org, novel_table in self.novel.items():
            novel_list = self.session.query(novel_table).all()
            org_keys = novel_table.c.keys()
            profile_list = self.session.query(self.profiles[org]).all()
            # Filter
            for novel in novel_list:
                filter_clauses = [
                    self.profiles[org].c[key] == getattr(novel, key)
                    for key in org_keys
                    if key not in ("ST", "clonal_complex", "species")
                ]
                exist = self.session.query(self.profiles[org]).filter(and_(*filter_clauses)).all()

                if exist:
                    exist = exist[0]
                    if sample == "":
                        onelap = prequery.filter(
                            and_(
                                Samples.ST == novel.ST,
                                Samples.organism == org,
                                Samples.ST <= -10,
                            )
                        ).all()
                    else:
                        onelap = prequery.filter(
                            and_(
                                Samples.ST == novel.ST,
                                Samples.organism == org,
                                Samples.ST <= -10,
                                Samples.CG_ID_sample == sample,
                            )
                        ).all()
                    for entry in onelap:
                        # review
                        if entry.pubmlst_ST == -1 and not overwrite:
                            self.logger.info(
                                f"Update: Sample {entry.CG_ID_sample} of organism {org}; Internal ST {novel.ST} is now linked to {exist.ST} '{exist}'"
                            )
                            self.update_sample(
                                {"CG_ID_sample": entry.CG_ID_sample},
                                {"pubmlst_ST": exist.ST},
                            )
                        # overwrite
                        elif overwrite:
                            self.logger.info(
                                f"Replacement: Sample {entry.CG_ID_sample} of organism {org}; Internal ST {novel.ST} is now {exist.ST} '{exist}'"
                            )
                            self.update_sample(
                                {"CG_ID_sample": entry.CG_ID_sample},
                                {"ST": exist.ST, "pubmlst_ST": exist.ST},
                            )

    def set_novel_ignored(self, sample=""):
        """Flags a sample as pubMLST resolved by merit of ignoring it"""
        query = self.session.query(Samples).filter(Samples.CG_ID_sample == sample).all()
        if len(query) > 0:
            self.logger.info(
                f"Ignore: Sample {query[0].CG_ID_sample} from organism {query[0].organism} with ST {query[0].ST}; is now flagged as resolved."
            )
            self.update_sample({"CG_ID_sample": query[0].CG_ID_sample}, {"pubmlst_ST": 0})
        else:
            self.logger.error(f"Sample {sample} not found in database. Verify name")

    def read_unresolved(self):
        """Lists all novel samples that current havent been flagged as resolved"""
        # ST currently not updated at all
        novelbkt = OrderedDict()
        prequery = (
            self.session.query(Samples)
            .filter(and_(Samples.ST <= -10, Samples.pubmlst_ST == -1))
            .all()
        )
        for entry in prequery:
            if entry.organism not in novelbkt:
                novelbkt[entry.organism] = dict()
            if entry.ST not in novelbkt[entry.organism]:
                novelbkt[entry.organism][entry.ST] = list()
            novelbkt[entry.organism][entry.ST].append(entry.CG_ID_sample)
        novelbkt = OrderedDict(sorted(novelbkt.items(), key=lambda t: t[0]))

        # ST updated on pubMLST but not marked as resolved:
        novelbkt2 = OrderedDict()
        postquery = (
            self.session.query(Samples)
            .filter(and_(Samples.ST <= -10, Samples.pubmlst_ST != -1, Samples.pubmlst_ST != 0))
            .all()
        )
        for entry in postquery:
            if entry.organism not in novelbkt2:
                novelbkt2[entry.organism] = dict()
            if entry.ST not in novelbkt2[entry.organism]:
                novelbkt2[entry.organism][entry.ST] = list()
            novelbkt2[entry.organism][entry.ST].append(entry.CG_ID_sample)

        # Unresolved samples and their respective error flags:
        novelbkt3 = OrderedDict()
        naquery = (
            self.session.query(Samples)
            .filter(and_(Samples.ST < 0, Samples.ST > -10, Samples.pubmlst_ST != 0))
            .all()
        )
        for entry in naquery:
            if entry.ST not in novelbkt3:
                novelbkt3[entry.ST] = dict()
            if entry.organism not in novelbkt3[entry.ST]:
                novelbkt3[entry.ST][entry.organism] = list()
            novelbkt3[entry.ST][entry.organism].append(entry.CG_ID_sample)
        novelbkt3 = OrderedDict(sorted(novelbkt3.items(), key=lambda t: t[0], reverse=True))

        codetrans = {
            -1: "Invalid pubMLST reference",
            -2: "Possibly novel allele, novel ST",
            -3: "Can't establish 7 loci due to low quality",
            -4: "Miscellaneous issues",
        }

        print("\n####Unresolved samples and their respective error flags:####\n")
        for k, v in novelbkt3.items():
            print(f"\n##Code {k} - {codetrans[k]}##")
            for x, y in v.items():
                if x is not None:
                    x = x.replace("_", " ").capitalize()
                print(f"{x} ({len(y)} samples):\n{sorted(y)}")
        if len(novelbkt3) == 0:
            print("None!")

        print("\n####ST updated on pubMLST but not marked as resolved:####\n")
        for k, v in novelbkt2.items():
            if k is not None:
                k = k.replace("_", " ").capitalize()
            print(f"Organism {k} ({len(v)}):")
            for x, y in v.items():
                print(f"{x}:{sorted(y)} ({len(y)} ST)")
        if len(novelbkt2) == 0:
            print("None!")

        print("\n####ST currently not updated at all:####\n")
        for k, v in novelbkt.items():
            if k is not None:
                k = k.replace("_", " ").capitalize()
            print(f"Organism {k} ({len(v)}):")
            for x, y in v.items():
                print(f"{x}:{sorted(y)} ({len(y)} novel ST)")
        if len(novelbkt) == 0:
            print("None!")

    def set_predictor(self, cg_sid: str, pks=dict()):
        """Helper function. Flags a set of seq_types as part of the final prediction.
        Uses optional pks[PK_NAME] = VALUE dictionary to distinguish in scenarios where an allele number has multiple hits
        """
        sample = self.session.query(Seq_types).filter(Seq_types.CG_ID_sample == cg_sid)

        if pks == dict():
            sample.update({Seq_types.st_predictor: 1})
        else:
            # Resets all previous predictors
            sample.update({Seq_types.st_predictor: None})
            # Set subset
            for loci, columns in pks.items():
                filter_clauses = [getattr(Seq_types, key) == val for key, val in columns.items()]
                sample.filter(and_(*filter_clauses)).update({Seq_types.st_predictor: 1})
        self.session.commit()

    def _build_allele_filter_clauses(self, alleles: dict, table) -> list:
        """Builds SQLAlchemy filter clauses for allele matching against a profile or novel table."""
        filter_clauses = []
        for key, val in alleles.items():
            col = table.c[key]
            if len(val) > 1:
                filter_clauses.append(or_(*[col == num for num in val]))
            else:
                filter_clauses.append(col == val[0])
        return filter_clauses

    def _query_st_profiles(self, alleles: dict, table) -> list:
        """Queries a profile or novel table with allele filter clauses and returns all matching rows."""
        filter_clauses = self._build_allele_filter_clauses(alleles, table)
        return self.session.query(table).filter(and_(*filter_clauses)).all()

    def _next_novel_st(self, organism: str) -> int:
        """Returns the next available negative novel ST: one below the current minimum, at most -10."""
        st = -9
        for entry in self.session.query(self.novel[organism]).all():
            if entry.ST < st:
                st = entry.ST
        return st - 1

    def _create_novel_st_entry(self, cg_sid: str, organism: str) -> int:
        """Creates a new novel ST row built from the sample's best alleles. Returns the new ST number."""
        st = self._next_novel_st(organism)
        best_alleles = self.read_best_alleles(cg_sid)
        new_entry: dict = {allele: columns["allele"] for allele, columns in best_alleles.items()}
        new_entry["ST"] = st
        self.add_rec(new_entry, self.novel[organism])
        return st

    def _allele_hit_score(self, allele) -> tuple:
        """Returns a comparable score tuple (span*identity, -evalue, contig_coverage) for one allele hit."""
        return (
            float(allele.span) * float(allele.identity),
            -float(allele.evalue),
            float(allele.contig_coverage),
        )

    def _score_profile(self, cg_sid: str, prof) -> tuple[dict, dict]:
        """For one profile row, fetches matching Seq_type alleles, keeps the best hit per locus,
        and returns (contig_names, score) where score sums spanid/eval/cc across all loci."""
        non_locus = {"ST", "clonal_complex", "species"}
        prof_keys = list(prof._fields)
        alleleconditions: list = []
        alleledict: dict = {}
        for index, allele_num in enumerate(prof):
            col_name = prof_keys[index]
            if col_name in non_locus:
                continue
            alleledict[col_name] = None
            alleleconditions.append(
                and_(Seq_types.loci == col_name, Seq_types.allele == allele_num)
            )
        all_alleles = (
            self.session.query(Seq_types)
            .filter(and_(Seq_types.CG_ID_sample == cg_sid, or_(*alleleconditions)))
            .all()
        )
        for allele in all_alleles:
            existing = alleledict[allele.loci]
            if existing is None or self._allele_hit_score(allele) > self._allele_hit_score(
                existing
            ):
                alleledict[allele.loci] = allele
        score: dict = {"spanid": 0.0, "eval": 0.0, "cc": 0.0}
        contig_names: dict = {}
        for locus, allele in alleledict.items():
            if allele is None:
                continue
            score["spanid"] += float(allele.span) * float(allele.identity)
            score["eval"] += float(allele.evalue)
            score["cc"] += float(allele.contig_coverage)
            contig_names[locus] = {"contig_name": str(allele.contig_name)}
        return contig_names, score

    def _pick_top_st(self, scores: dict) -> int | str:
        """Selects the ST with the highest composite score (spanid → eval → contig_coverage)."""
        top_st: int | str | None = None
        top_spanid = -1.0
        top_eval = float("inf")
        top_cc = -1.0
        for st, val in scores.items():
            if (
                val["spanid"] > top_spanid
                or (val["spanid"] == top_spanid and val["eval"] < top_eval)
                or (val["spanid"] == top_spanid and val["eval"] == top_eval and val["cc"] > top_cc)
            ):
                top_spanid = val["spanid"]
                top_eval = val["eval"]
                top_cc = val["cc"]
                top_st = st
        return top_st

    def read_st(self, cg_sid: str) -> int | str:
        """Takes a CG_ID_sample and predicts the correct ST."""
        organism: str | None = (
            self.session.query(Samples.organism).filter(Samples.CG_ID_sample == cg_sid).scalar()
        )
        if organism is None:
            self.logger.warning(
                f"No organism set for {cg_sid}. Most likely control sample. Setting ST to -1"
            )
            return -1

        threshold = True
        alleles, allelediff = self.read_unique_alleles(cg_sid, organism, threshold)
        if allelediff < 0:
            threshold = False
            alleles, allelediff = self.read_unique_alleles(cg_sid, organism, threshold)
            if allelediff < 0:
                self.logger.warning(
                    f"Insufficient allele hits to establish ST for sample {cg_sid}, even without thresholds. Setting ST to -3"
                )
                self.set_predictor(cg_sid)
                return -3

        # Try matching against the curated profile table
        output = self._query_st_profiles(alleles, self.profiles[organism])
        if output:
            st_list = [row.ST for row in output]
            best = self.read_best_st(cg_sid, st_list, "profile")
            if threshold and len(st_list) > 1:
                self.logger.warning(
                    f"Multiple ST within threshold found for sample {cg_sid}, list: {st_list}. Established ST{best} as best hit."
                )
            return best

        # Try matching against the novel ST table (only when hits are above threshold)
        if threshold:
            self.logger.info(
                f"Sample {cg_sid} on {organism} has novel ST reliably established. Searching for prior novel definition..."
            )
            output = self._query_st_profiles(alleles, self.novel[organism])
            if output:
                st_list = [row.ST for row in output]
                best = self.read_best_st(cg_sid, st_list, "novel")
                if len(st_list) > 1:
                    self.logger.warning(
                        f"Multiple ST within novel threshold found for sample {cg_sid}, list: {st_list}. Established ST{best} as best hit."
                    )
                return best
            # No prior novel match — create a new novel ST
            new_st = self._create_novel_st_entry(cg_sid, organism)
            return self.read_best_st(cg_sid, [new_st], "novel")

        self.logger.warning(
            f"Sample {cg_sid} on {organism} has an allele set but hits are low-quality and do not resolve to an ST. Setting ST to -2"
        )
        best_set = self.read_best_alleles(cg_sid)
        self.set_predictor(cg_sid, best_set)
        return -2

    def read_best_st(self, cg_sid: str, st_list: list, type: str = "profile") -> int | str:
        """Establishes which ST is most likely by criteria id*span → eval → contig coverage
        and flags the involved alleles as predictors."""
        organism: str = (
            self.session.query(Samples.organism).filter(Samples.CG_ID_sample == cg_sid).scalar()
        )
        table: Table = self.profiles[organism] if type == "profile" else self.novel[organism]
        profiles = [self.session.query(table).filter(text(f"ST={st}")).first() for st in st_list]
        scores: dict = {}
        best_alleles: dict = {}
        for prof in profiles:
            contig_names, score = self._score_profile(cg_sid, prof)
            scores[prof.ST] = score
            best_alleles[prof.ST] = contig_names
        top_st = self._pick_top_st(scores)
        self.set_predictor(cg_sid, best_alleles[top_st])
        return top_st

    def read_best_alleles(self, cg_sid: str):
        """Establishes which allele set (for bad samples) is most likely by criteria span* id -> eval -> contig coverage"""
        hits: list[tuple[str, str, float, float, float, float, str]] = (
            self.session.query(
                Seq_types.contig_name,
                Seq_types.loci,
                Seq_types.span,
                Seq_types.identity,
                Seq_types.evalue,
                Seq_types.contig_coverage,
                Seq_types.allele,
            )
            .filter(Seq_types.CG_ID_sample == cg_sid)
            .all()
        )
        bestHits: dict[str, dict[str, str]] = {}
        alleledict: dict[str, list[float]] = {}
        for allele in hits:
            if allele.loci in bestHits:
                if (
                    (
                        allele.identity * allele.span
                        > alleledict[allele.loci][0] * alleledict[allele.loci][3]
                    )
                    or (
                        allele.identity * allele.span
                        == alleledict[allele.loci][0] * alleledict[allele.loci][3]
                        and float(allele.evalue) < float(alleledict[allele.loci][1])
                    )
                    or (
                        allele.identity * allele.span
                        == alleledict[allele.loci][0] * alleledict[allele.loci][3]
                        and float(allele.evalue) == float(alleledict[allele.loci][1])
                        and allele.contig_coverage > alleledict[allele.loci][2]
                    )
                ):
                    bestHits[allele.loci]["contig_name"] = allele.contig_name
                    alleledict[allele.loci] = [
                        allele.identity,
                        allele.evalue,
                        allele.contig_coverage,
                        allele.span,
                    ]
            else:
                bestHits[allele.loci] = {
                    "contig_name": allele.contig_name,
                    "allele": allele.allele,
                }
                alleledict[allele.loci] = [
                    allele.identity,
                    allele.evalue,
                    allele.contig_coverage,
                    allele.span,
                ]
        return bestHits

    def read_unique_alleles(self, cg_sid: str, organism: str, threshold=True):
        """Returns a dict containing all unique alleles at every loci, and allele difference from expected"""
        tspan = (self.threshold.mlst_span) / 100.0
        if threshold:
            tid = float(self.threshold.mlst_id)
            hits: list[tuple[str, str]] = (
                self.session.query(Seq_types.loci, Seq_types.allele)
                .filter(
                    Seq_types.CG_ID_sample == cg_sid,
                    Seq_types.identity >= tid,
                    Seq_types.span >= tspan,
                )
                .all()
            )
        else:
            hits: list[tuple[str, str]] = (
                self.session.query(Seq_types.loci, Seq_types.allele)
                .filter(Seq_types.CG_ID_sample == cg_sid)
                .all()
            )

        # Establish number of unique hits
        uniqueDict: dict[str, list[str]] = {}
        for hit in hits:
            if hit.loci not in uniqueDict.keys():
                uniqueDict[hit.loci] = [hit.allele]
            elif hit.allele not in uniqueDict[hit.loci]:
                uniqueDict[hit.loci].append(hit.allele)
        non_allele_columns = 1
        if "clonal_complex" in self.profiles[organism].columns.keys():
            non_allele_columns += 1
        if "species" in self.profiles[organism].columns.keys():
            non_allele_columns += 1
        allele_overabundance = len(uniqueDict.keys()) - (
            len(self.profiles[organism].columns.values()) - non_allele_columns
        )
        return [uniqueDict, allele_overabundance]
