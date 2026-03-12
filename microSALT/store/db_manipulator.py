"""Delivers and fetches data from the database
By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import hashlib
import sys
import warnings
from collections import OrderedDict
from datetime import datetime, timezone

from dateutil.parser import parse
from sqlalchemy import inspect as sa_inspect, MetaData, desc, or_, and_, text

from microSALT import __version__
from microSALT.config import Folders, Threshold
from microSALT.exc.exceptions import RefUpdateLockError
from microSALT.store.database import get_session, get_engine
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
                self.add_rec(
                    {"name": f"profile_{k}", "version": "0"},
                    "Versions",
                    force=True,
                )
                self.logger.info(f"Profile table profile_{k} created and populated")
        for k, v in self.novel.items():
            if not inspector.has_table(f"novel_{k}"):
                self.novel[k].create(self.engine)
                self.add_rec(
                    {"name": f"novel_{k}", "version": "0"},
                    "Versions",
                    force=True,
                )
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

    def add_rec(self, data_dict: dict[str, str], tablename: str, force=False):
        """Adds a record to the specified table through a dict with columns as keys."""
        pk_list = list()
        # Non-orm
        if not isinstance(tablename, str):
            # check for existence
            table = tablename
            pk_list = table.primary_key.columns.keys()
            filter_clauses = [table.c[pk] == data_dict[pk] for pk in pk_list]
            exist = self.session.query(table).filter(or_(*filter_clauses)).all()
            # Add record
            if len(exist) == 0:
                data = table.insert()
                # Loads any dates as datetime objects
                for k, v in data_dict.items():
                    if isinstance(v, str):
                        try:
                            parse(v, fuzzy=False)
                            data_dict[k] = datetime.strptime(v, "%Y-%m-%d %H:%M:%S")
                        except ValueError as ve:
                            if len(ve.args) > 0 and ve.args[0].startswith(
                                "unconverted data remains: "
                            ):
                                data_dict[k] = datetime.strptime(v, "%Y-%m-%d %H:%M:%S.%f")
                            else:
                                pass
                self.session.execute(data, data_dict)
                self.session.commit()
                self.logger.info(f"Added entry to table {tablename.fullname}")
        # ORM
        else:
            try:
                table = _resolve_orm_table(tablename)
                # Check for existing entry
                pk_list = table.__table__.primary_key.columns.keys()
            except KeyError:
                self.logger.error(
                    f"Attempted to access table {tablename} which has not been created"
                )
                return
            pk_values = list()
            for item in pk_list:
                pk_values.append(data_dict[item])
            existing = self.session.get(table, pk_values)
            # Add record
            if not existing or force:
                newobj = table()
                # Loads any dates as datetime objects
                for k, v in data_dict.items():
                    if isinstance(v, str):
                        try:
                            parse(v, fuzzy=False)
                            data_dict[k] = datetime.strptime(v, "%Y-%m-%d %H:%M:%S")
                        except ValueError as ve:
                            if len(ve.args) > 0 and ve.args[0].startswith(
                                "unconverted data remains: "
                            ):
                                data_dict[k] = datetime.strptime(v, "%Y-%m-%d %H:%M:%S.%f")
                            else:
                                pass
                for k, v in data_dict.items():
                    setattr(newobj, k, v)
                self.session.add(newobj)
                self.session.commit()
            else:
                self.logger.warning(
                    f"Record [{', '.join(pk_list)}]=[{', '.join(pk_values)}] in table {tablename} already exists"
                )

    def upd_rec(self, req_dict: dict[str, str], tablename: str, upd_dict: dict[str, str]):
        """Updates a record to the specified table through a dict with columns as keys."""
        table = _resolve_orm_table(tablename)
        self.logger.debug(f"Updating table {tablename} with {upd_dict}")
        filter_clauses = [getattr(table, k) == v for k, v in req_dict.items() if v is not None]
        query = self.session.query(table).filter(and_(*filter_clauses))
        if len(query.all()) > 1:
            self.logger.error("More than 1 record found when orm updating. Exited.")
            sys.exit()
        else:
            # If the primary key CG_ID_sample is being renamed, propagate the
            # change to child tables first (bulk query.update bypasses ORM
            # cascade logic and would otherwise leave orphaned rows).
            if tablename == "Samples" and "CG_ID_sample" in upd_dict:
                old_id = req_dict.get("CG_ID_sample")
                new_id = upd_dict["CG_ID_sample"]
                if old_id and old_id != new_id:
                    for child_table in (Seq_types, Resistances, Expacs, Collections):
                        self.session.query(child_table).filter(
                            child_table.CG_ID_sample == old_id
                        ).update({"CG_ID_sample": new_id})
            query.update(upd_dict)
            self.session.commit()
        self.logger.debug(f"Updated table {tablename} with {upd_dict} for {req_dict}")

    def purge_rec(self, name: str, type: str):
        """Removes seq_data, resistances, sample(s) and possibly project"""
        entries = list()
        if type == "Projects":
            entries.append(
                self.session.query(Expacs).filter(Expacs.CG_ID_sample.like(f"{name}%")).all()
            )
            entries.append(
                self.session.query(Seq_types).filter(Seq_types.CG_ID_sample.like(f"{name}%")).all()
            )
            entries.append(
                self.session.query(Resistances)
                .filter(Resistances.CG_ID_sample.like(f"{name}%"))
                .all()
            )
            entries.append(
                self.session.query(Samples).filter(Samples.CG_ID_sample.like(f"{name}%")).all()
            )
            # entries.append(self.session.query(Projects).filter(Projects.CG_ID_project==name).all())
        elif type == "Samples":
            entries.append(self.session.query(Expacs).filter(Expacs.CG_ID_sample == name).all())
            entries.append(
                self.session.query(Seq_types).filter(Seq_types.CG_ID_sample == name).all()
            )
            entries.append(
                self.session.query(Resistances).filter(Resistances.CG_ID_sample == name).all()
            )
            entries.append(self.session.query(Samples).filter(Samples.CG_ID_sample == name).all())
        elif type == "Collections":
            entries.append(
                self.session.query(Collections).filter(Collections.ID_collection == name).all()
            )
        else:
            self.logger.error(f"Incorrect type {type} specified for removal of {name}. Check code")
            sys.exit()
        for entry in entries:
            for instance in entry:
                self.session.delete(instance)
        self.session.commit()
        self.logger.info(f"Removed information for {name}")

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

    def read_exists(self, table: str, item: dict[str, str]):
        """Takes a k-v pair and checks for the entrys existence in the given table"""
        orm_table = _resolve_orm_table(table)
        filter_clauses = [getattr(orm_table, k) == v for k, v in item.items()]
        entry = self.session.query(orm_table).filter(and_(*filter_clauses)).scalar()
        return entry is not None

    def read_version(self, name: str):
        """Gets the version from a given name. Should be generalized to return any value for any input"""
        version = self.session.query(Versions).filter(Versions.name == name).scalar()
        if version is None:
            return "0"
        else:
            return version.version

    def read_report(self, name: str):
        # Sort based on version
        prev_report = []
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
        totalstring = list()
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

        prev_report = self.read_report(name)
        # Compare
        if prev_report:
            if "steps_aggregate" in dir(prev_report) and prev_report.steps_aggregate != hashstring:
                self.add_rec(
                    {
                        "CG_ID_project": name,
                        "steps_aggregate": hashstring,
                        "date": dt,
                        "version": prev_report.version + 1,
                    },
                    "Reports",
                )
        else:
            self.add_rec(
                {
                    "CG_ID_project": name,
                    "steps_aggregate": hashstring,
                    "date": dt,
                    "version": 1,
                },
                "Reports",
            )

    def sync_novel(self, overwrite=False, sample=""):
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
                            self.upd_rec(
                                {"CG_ID_sample": entry.CG_ID_sample},
                                "Samples",
                                {"pubmlst_ST": exist.ST},
                            )
                        # overwrite
                        elif overwrite:
                            self.logger.info(
                                f"Replacement: Sample {entry.CG_ID_sample} of organism {org}; Internal ST {novel.ST} is now {exist.ST} '{exist}'"
                            )
                            self.upd_rec(
                                {"CG_ID_sample": entry.CG_ID_sample},
                                "Samples",
                                {"ST": exist.ST, "pubmlst_ST": exist.ST},
                            )

    def rm_novel(self, sample=""):
        """Flags a sample as pubMLST resolved by merit of ignoring it"""
        query = self.session.query(Samples).filter(Samples.CG_ID_sample == sample).all()
        if len(query) > 0:
            self.logger.info(
                f"Ignore: Sample {query[0].CG_ID_sample} from organism {query[0].organism} with ST {query[0].ST}; is now flagged as resolved."
            )
            self.upd_rec({"CG_ID_sample": query[0].CG_ID_sample}, "Samples", {"pubmlst_ST": 0})
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

    def setPredictor(self, cg_sid: str, pks=dict()):
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

    def read_st(self, cg_sid: str):
        """Takes a CG_ID_sample and predicts the correct ST"""
        threshold = True
        organism = (
            self.session.query(Samples.organism).filter(Samples.CG_ID_sample == cg_sid).scalar()
        )
        if organism is None:
            self.logger.warning(
                f"No organism set for {cg_sid}. Most likely control sample. Setting ST to -1"
            )
            return -1
        [alleles, allelediff] = self.read_unique_alleles(cg_sid, organism, threshold)
        if allelediff < 0:
            threshold = False
            [alleles, allelediff] = self.read_unique_alleles(cg_sid, organism, threshold)
            if allelediff < 0:
                self.logger.warning(
                    f"Insufficient allele hits to establish ST for sample {cg_sid}, even without thresholds. Setting ST to -3"
                )
                self.setPredictor(cg_sid)
                return -3

        # Tests all allele combinations found to see if any of them result in ST
        filter_clauses = []
        for key, val in alleles.items():
            col = self.profiles[organism].c[key]
            if len(val) > 1:
                filter_clauses.append(or_(*[col == num for num in val]))
            else:
                filter_clauses.append(col == val[0])
        output = self.session.query(self.profiles[organism]).filter(and_(*filter_clauses)).all()

        # Check for existence in profile database
        if len(output) > 1:
            STlist = list()
            for st in output:
                STlist.append(st.ST)
            best = self.read_best_st(cg_sid, STlist, "profile")
            if threshold:
                self.logger.warning(
                    f"Multiple ST within threshold found for sample {cg_sid}, list: {STlist}. Established ST{best} as best hit."
                )
            return best
        elif len(output) == 1:
            # Arbitary call
            return self.read_best_st(cg_sid, [output[0].ST], "profile")
        # Check for existence in novel database
        elif threshold:
            self.logger.info(
                f"Sample {cg_sid} on {organism} has novel ST reliably established. Searching for prior novel definition..."
            )
            filter_clauses = []
            for key, val in alleles.items():
                col = self.novel[organism].c[key]
                if len(val) > 1:
                    filter_clauses.append(or_(*[col == num for num in val]))
                else:
                    filter_clauses.append(col == val[0])
            output = self.session.query(self.novel[organism]).filter(and_(*filter_clauses)).all()

            if len(output) > 1:
                STlist = list()
                for st in output:
                    STlist.append(st.ST)
                best = self.read_best_st(cg_sid, STlist, "novel")
                if threshold:
                    self.logger.warning(
                        f"Multiple ST within novel threshold found for sample {cg_sid}, list: {STlist}. Established ST{best} as best hit."
                    )
                return best
            elif len(output) == 1:
                return self.read_best_st(cg_sid, [output[0].ST], "novel")
            else:
                # Create new novel ST
                # Set ST -10 per default, or one below the current min, whichever is smaller.
                st = -9
                query = self.session.query(self.novel[organism]).all()
                for entry in query:
                    if entry.ST < st:
                        st = entry.ST
                st = st - 1

                bestSet = self.read_best_alleles(cg_sid)
                newEntry = dict()
                for allele, columns in bestSet.items():
                    newEntry[allele] = columns["allele"]
                newEntry["ST"] = st
                self.add_rec(newEntry, self.novel[organism])
                return self.read_best_st(cg_sid, [st], "novel")
        else:
            self.logger.warning(
                f"Sample {cg_sid} on {organism} has an allele set but hits are low-quality and do not resolve to an ST. Setting ST to -2"
            )
            bestSet = self.read_best_alleles(cg_sid)
            self.setPredictor(cg_sid, bestSet)
            return -2

    def read_best_st(self, cg_sid: str, st_list: list, type="profile"):
        """Takes in a list of ST and a sample.
        Establishes which ST is most likely by criteria id*span -> eval -> contig coverage
        & flags involved alleles"""
        profiles = list()
        scores = dict()
        bestalleles = dict()
        organism = (
            self.session.query(Samples.organism).filter(Samples.CG_ID_sample == cg_sid).scalar()
        )
        for st in st_list:
            scores[st] = dict()
            bestalleles[st] = dict()
            scores[st]["spanid"] = 0
            scores[st]["eval"] = 0
            scores[st]["cc"] = 0
            scores[st]["span"] = 0
            if type == "profile":
                profiles.append(
                    self.session.query(self.profiles[organism]).filter(text(f"ST={st}")).first()
                )
            elif type == "novel":
                profiles.append(
                    self.session.query(self.novel[organism]).filter(text(f"ST={st}")).first()
                )

        # Get values for each allele set that resolves an ST
        for prof in profiles:
            prof_keys = list(prof._fields)
            alleleconditions = list()
            alleledict = dict()

            for index, allele in enumerate(prof):
                if (
                    "ST" not in prof_keys[index]
                    and "clonal_complex" not in prof_keys[index]
                    and "species" not in prof_keys[index]
                ):
                    alleledict[prof_keys[index]] = ""
                    alleleconditions.append(
                        and_(
                            Seq_types.loci == prof_keys[index],
                            Seq_types.allele == allele,
                        )
                    )

            all_alleles = (
                self.session.query(Seq_types)
                .filter(and_(Seq_types.CG_ID_sample == cg_sid, or_(*alleleconditions)))
                .all()
            )

            # Keep only best hit each loci
            for allele in all_alleles:
                if alleledict[allele.loci] == "":
                    alleledict[allele.loci] = allele
                else:
                    old_al = alleledict[allele.loci]

                    if allele.span * allele.identity >= old_al.span * old_al.identity:
                        if allele.span * allele.identity > old_al.span * old_al.identity:
                            alleledict[allele.loci] = allele
                        elif float(allele.evalue) <= float(old_al.evalue):
                            if float(allele.evalue) < float(old_al.evalue):
                                alleledict[allele.loci] = allele
                            elif allele.contig_coverage > old_al.contig_coverage:
                                alleledict[allele.loci] = allele

            # Create score dict for the ST
            for key, allele in alleledict.items():
                scores[prof.ST]["spanid"] += allele.span * allele.identity
                scores[prof.ST]["eval"] += float(allele.evalue)
                scores[prof.ST]["cc"] += allele.contig_coverage
                if allele.loci not in bestalleles[prof.ST].keys():
                    bestalleles[prof.ST][allele.loci] = dict()
                if "contig_name" not in bestalleles[prof.ST][allele.loci].keys():
                    bestalleles[prof.ST][allele.loci]["contig_name"] = str(allele.contig_name)

        # Establish best ST
        topST = ""
        topID = 0
        topEval = 100
        topCC = 0
        for key, val in scores.items():
            if scores[key]["spanid"] > topID:
                topID = scores[key]["spanid"]
                topEval = scores[key]["eval"]
                topCC = scores[key]["cc"]
                topST = key
            elif scores[key]["spanid"] == topID and scores[key]["eval"] < topEval:
                topID = scores[key]["spanid"]
                topEval = scores[key]["eval"]
                topCC = scores[key]["cc"]
                topST = key
            elif (
                scores[key]["spanid"] == topID
                and scores[key]["eval"] == topEval
                and scores[key]["cc"] > topCC
            ):
                topID = scores[key]["spanid"]
                topEval = scores[key]["eval"]
                topCC = scores[key]["cc"]
                topST = key
        self.setPredictor(cg_sid, bestalleles[topST])
        return topST

    def read_best_alleles(self, cg_sid: str):
        """Establishes which allele set (for bad samples) is most likely by criteria span* id -> eval -> contig coverage"""
        hits = (
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
        bestHits = dict()
        alleledict = dict()
        for allele in hits:
            if allele.loci not in bestHits.keys():
                bestHits[allele.loci] = dict()
                bestHits[allele.loci]["contig_name"] = allele.contig_name
                bestHits[allele.loci]["allele"] = allele.allele
                alleledict[allele.loci] = [
                    allele.identity,
                    allele.evalue,
                    allele.contig_coverage,
                    allele.span,
                ]
            else:
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
        return bestHits

    def read_unique_alleles(self, cg_sid: str, organism: str, threshold=True):
        """Returns a dict containing all unique alleles at every loci, and allele difference from expected"""
        tid = float(self.threshold.mlst_id)
        tspan = (self.threshold.mlst_span) / 100.0
        if threshold:
            hits = (
                self.session.query(Seq_types.loci, Seq_types.allele)
                .filter(
                    Seq_types.CG_ID_sample == cg_sid,
                    Seq_types.identity >= tid,
                    Seq_types.span >= tspan,
                )
                .all()
            )
        else:
            hits = (
                self.session.query(Seq_types.loci, Seq_types.allele)
                .filter(Seq_types.CG_ID_sample == cg_sid)
                .all()
            )

        # Establish number of unique hits
        uniqueDict = dict()
        for hit in hits:
            if hit.loci not in uniqueDict.keys():
                uniqueDict[hit.loci] = list()
                uniqueDict[hit.loci].append(hit.allele)
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
