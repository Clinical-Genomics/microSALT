""" Delivers and fetches data from the database
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import hashlib
import sys
import warnings

from collections import OrderedDict
from datetime import datetime, timezone

from dateutil.parser import parse
from sqlalchemy import inspect as sa_inspect, MetaData, desc, create_engine, or_, and_, text
from sqlalchemy.orm import sessionmaker

# maintain the same connection per thread
from sqlalchemy.pool import SingletonThreadPool

from microSALT import __version__
from microSALT.store.orm_models import (
    app,
    Collections,
    Expacs,
    Projects,
    Reports,
    Resistances,
    Samples,
    Seq_types,
    Versions,
)
from microSALT.store.models import Profiles, Novel


class DB_Manipulator:
    def __init__(self, config, log):
        self.config = config
        self.logger = log
        self.engine = create_engine(
            app.config["SQLALCHEMY_DATABASE_URI"], poolclass=SingletonThreadPool
        )
        Session = sessionmaker(bind=self.engine)
        self.session = Session()
        self.metadata = MetaData(self.engine)
        self.profiles = Profiles(self.metadata, self.config, self.logger).tables
        self.novel = Novel(self.metadata, self.config, self.logger).tables
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
        for k, v in self.profiles.items():
            if not inspector.has_table(f"profile_{k}"):
                self.profiles[k].create()
                self.init_profiletable(k, v)
                self.add_rec(
                    {"name": f"profile_{k}", "version": "0"},
                    "Versions",
                    force=True,
                )
                self.logger.info(f"Profile table profile_{k} initialized")
        for k, v in self.novel.items():
            if not inspector.has_table(f"novel_{k}"):
                self.novel[k].create()
                self.add_rec(
                    {"name": f"novel_{k}", "version": "0"},
                    "Versions",
                    force=True,
                )
                self.logger.info(f"Profile table novel_{k} initialized")

    def add_rec(self, data_dict: dict[str, str], tablename: str, force=False):
        """Adds a record to the specified table through a dict with columns as keys."""
        pk_list = list()
        # Non-orm
        if not isinstance(tablename, str):
            # check for existence
            table = tablename
            pk_list = table.primary_key.columns.keys()
            args = list()
            for pk in pk_list:
                args.append(f"table.c.{pk}=={data_dict[pk]}")
            args = "or_(" + ",".join(args) + ")"
            exist = self.session.query(table).filter(eval(args)).all()
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
                with self.engine.connect() as conn:
                    conn.execute(data, data_dict)
                    conn.commit()
                self.logger.info(f"Added entry to table {tablename.fullname}")
        # ORM
        else:
            try:
                table = eval(tablename)
                # Check for existing entry
                pk_list = table.__table__.primary_key.columns.keys()
            except Exception as e:
                self.logger.error(
                    f"Attempted to access table {tablename} which has not been created"
                )
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
        table = eval(tablename)
        self.logger.debug(f"Updating table {tablename} with {upd_dict}")
        argy = list()
        for k, v in req_dict.items():
            if v != None:
                argy.append(f".filter(table.{k}=='{v}')")
        filter = "".join(argy)
        megastring = f"self.session.query(table){filter}"
        if len(eval(megastring + ".all()")) > 1:
            self.logger.error("More than 1 record found when orm updating. Exited.")
            sys.exit()
        else:
            eval(megastring + ".update(upd_dict)")
            self.session.commit()
        self.logger.debug(f"Updated table {tablename} with {upd_dict} for {req_dict}")

    def purge_rec(self, name: str, type: str):
        """Removes seq_data, resistances, sample(s) and possibly project"""
        entries = list()
        if type == "Projects":
            entries.append(
                self.session.query(Expacs)
                .filter(Expacs.CG_ID_sample.like(f"{name}%"))
                .all()
            )
            entries.append(
                self.session.query(Seq_types)
                .filter(Seq_types.CG_ID_sample.like(f"{name}%"))
                .all()
            )
            entries.append(
                self.session.query(Resistances)
                .filter(Resistances.CG_ID_sample.like(f"{name}%"))
                .all()
            )
            entries.append(
                self.session.query(Samples)
                .filter(Samples.CG_ID_sample.like(f"{name}%"))
                .all()
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
            self.logger.error(
                f"Incorrect type {type} specified for removal of {name}. Check code"
            )
            sys.exit()
        for entry in entries:
            for instance in entry:
                self.session.delete(instance)
                self.session.commit()
        self.logger.info(f"Removed information for {name}")

    def query_rec(self, tablename: str, filters: dict[str, str]):
        """Fetches records table, using a primary-key dict with columns as keys.
        Non-PK are ignored"""
        # Non-orm
        if not isinstance(tablename, str):
            # check for existence
            table = tablename
            pk_list = table.primary_key.columns.keys()
            args = list()
            for k, v in filters.items():
                args.append(f"table.c.{k}=={v}")
            args = "or_(" + ",".join(args) + ")"
            exist = self.session.query(table).filter(eval(args)).all()
            return exist
        # ORM
        else:
            table = eval(tablename)
            args = list()
            for k, v in filters.items():
                if v != None:
                    args.append(f"table.{k}=='{v}'")
            filter = " and ".join(args)
            entries = self.session.query(table).filter(eval(filter)).all()
            return entries

    def top_index(self, table_str: str, filters: dict[str, str], column: str):
        """Fetches the top index from column of table, by applying a dict with columns as keys."""
        table = eval(table_str)
        args = list()
        for k, v in filters.items():
            if v != None:
                args.append(f"table.{k}=='{v}'")
        filter = " and ".join(args)
        entry = (
            self.session.query(table)
            .filter(eval(filter))
            .order_by(desc(eval(f"{table_str}.{column}")))
            .limit(1)
            .all()
        )
        if entry == []:
            return int(-1)
        else:
            return eval(f"entry[0].{column}")

    def reload_profiletable(self, organism: str):
        """Drop the named non-orm table, then load it with fresh data"""
        table = self.profiles[organism]
        self.logger.debug(f"Reloading profile table for {organism}")
        self.profiles[organism].drop()
        self.logger.debug(f"Dropped profile table for {organism}")
        self.profiles[organism].create()
        self.logger.debug(f"Recreated profile table for {organism}")
        self.init_profiletable(organism, table)
        self.logger.debug(f"Initialized profile table for {organism}")

    def init_profiletable(self, filename: str, table):
        """Creates profile tables by looping, since a lot of infiles exist"""
        data = table.insert()
        linedict = dict.fromkeys(table.c.keys())
        file_path = f"{self.config['folders']['profiles']}/{filename}"
        self.logger.debug(f"Opening profile file: {file_path}")
        with open(file_path, "r") as fh:
            # Skips header
            head = fh.readline()
            head = head.rstrip().split("\t")
            self.logger.debug(f"Header columns: {head}")
            line_num = 0
            for line in fh:
                line_num += 1
                line = line.rstrip().split("\t")
                self.logger.debug(f"Processing line {line_num}: {line}")
                index = 0
                while index < len(line):
                    linedict[head[index]] = line[index]
                    index = index + 1
                self.logger.debug(f"Linedict before insert: {linedict}")
                try:
                    with self.engine.connect() as conn:
                        conn.execute(data, linedict)
                        conn.commit()
                    self.logger.debug(f"Inserted line {line_num} into table")
                except Exception as e:
                    self.logger.error(f"Failed to insert line {line_num}: {e}")
        self.logger.debug(f"Initialized profile table for {filename} with {line_num} entries")

    def get_columns(self, tablename: str):
        """Returns all records for a given ORM table"""
        table = eval(tablename)
        return dict.fromkeys(table.__table__.columns.keys())

    def exists(self, table, item: dict[str, str]):
        """Takes a k-v pair and checks for the entrys existence in the given table"""
        filterstring = ""
        for k, v in item.items():
            filterstring += f"{table}.{k}=='{v}',"
        filterstring = filterstring[:-1]
        table = eval(table)
        entry = self.session.query(table).filter(eval(filterstring)).scalar()
        if entry is None:
            return False
        else:
            return True

    def get_version(self, name: str):
        """Gets the version from a given name. Should be generalized to return any value for any input"""
        version = self.session.query(Versions).filter(Versions.name == name).scalar()
        if version is None:
            return "0"
        else:
            return version.version

    def get_report(self, name: str):
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

        prev_report = self.get_report(name)
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
                args = list()
                for key in org_keys:
                    if key != "ST" and key != "clonal_complex" and key != "species":
                        args.append(
                            f"self.profiles[org].c.{key}=={eval(f'novel.{key}')}"
                        )
                args = "and_(" + ",".join(args) + ")"
                exist = self.session.query(self.profiles[org]).filter(eval(args)).all()

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

    def list_unresolved(self):
        """Lists all novel samples that current havent been flagged as resolved"""
        # ST currently not updated at all
        novelbkt = OrderedDict()
        prequery = (
            self.session.query(Samples)
            .filter(and_(Samples.ST <= -10, Samples.pubmlst_ST == -1))
            .all()
        )
        for entry in prequery:
            if not entry.organism in novelbkt:
                novelbkt[entry.organism] = dict()
            if not entry.ST in novelbkt[entry.organism]:
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
            if not entry.organism in novelbkt2:
                novelbkt2[entry.organism] = dict()
            if not entry.ST in novelbkt2[entry.organism]:
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
            if not entry.ST in novelbkt3:
                novelbkt3[entry.ST] = dict()
            if not entry.organism in novelbkt3[entry.ST]:
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
                print(
                    f"{x} ({len(y)} samples):\n{sorted(y)}"
                )
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
                arglist = list()
                for key, val in columns.items():
                    arglist.append(f"Seq_types.{key}=='{val}'")
                    args = "and_(" + ", ".join(arglist) + ")"
                sample.filter(eval(args)).update({Seq_types.st_predictor: 1})
        self.session.commit()

    def alleles2st(self, cg_sid: str):
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
        [alleles, allelediff] = self.get_unique_alleles(cg_sid, organism, threshold)
        if allelediff < 0:
            threshold = False
            [alleles, allelediff] = self.get_unique_alleles(cg_sid, organism, threshold)
            if allelediff < 0:
                self.logger.warning(
                    f"Insufficient allele hits to establish ST for sample {cg_sid}, even without thresholds. Setting ST to -3"
                )
                self.setPredictor(cg_sid)
                return -3

        # Tests all allele combinations found to see if any of them result in ST
        filter = list()
        for key, val in alleles.items():
            subfilter = list()
            for num in val:
                subfilter.append(f" self.profiles[organism].c.{key}=={num} ")
            subfilter = ",".join(subfilter)
            if len(val) > 1:
                subfilter = f"or_({subfilter})"
            filter.append(subfilter)
        filter = ",".join(filter)
        filter = f"and_({filter})"
        output = self.session.query(self.profiles[organism]).filter(eval(filter)).all()

        # Check for existence in profile database
        if len(output) > 1:
            STlist = list()
            for st in output:
                STlist.append(st.ST)
            best = self.bestST(cg_sid, STlist, "profile")
            if threshold:
                self.logger.warning(
                    f"Multiple ST within threshold found for sample {cg_sid}, list: {STlist}. Established ST{best} as best hit."
                )
            return best
        elif len(output) == 1:
            # Arbitary call
            return self.bestST(cg_sid, [output[0].ST], "profile")
        # Check for existence in novel database
        elif threshold:
            self.logger.info(
                f"Sample {cg_sid} on {organism} has novel ST reliably established. Searching for prior novel definition..."
            )
            filter = list()
            for key, val in alleles.items():
                subfilter = list()
                for num in val:
                    subfilter.append(f" self.novel[organism].c.{key}=={num} ")
                subfilter = ",".join(subfilter)
                if len(val) > 1:
                    subfilter = f"or_({subfilter})"
                filter.append(subfilter)
            filter = ",".join(filter)
            filter = f"and_({filter})"
            output = self.session.query(self.novel[organism]).filter(eval(filter)).all()

            if len(output) > 1:
                STlist = list()
                for st in output:
                    STlist.append(st.ST)
                best = self.bestST(cg_sid, STlist, "novel")
                if threshold:
                    self.logger.warning(
                        f"Multiple ST within novel threshold found for sample {cg_sid}, list: {STlist}. Established ST{best} as best hit."
                    )
                return best
            elif len(output) == 1:
                return self.bestST(cg_sid, [output[0].ST], "novel")
            else:
                # Create new novel ST
                # Set ST -10 per default, or one below the current min, whichever is smaller.
                st = -9
                query = self.session.query(self.novel[organism]).all()
                for entry in query:
                    if entry.ST < st:
                        st = entry.ST
                st = st - 1

                bestSet = self.bestAlleles(cg_sid)
                newEntry = dict()
                for allele, columns in bestSet.items():
                    newEntry[allele] = columns["allele"]
                newEntry["ST"] = st
                self.add_rec(newEntry, self.novel[organism])
                return self.bestST(cg_sid, [st], "novel")
        else:
            self.logger.warning(
                f"Sample {cg_sid} on {organism} has an allele set but hits are low-quality and do not resolve to an ST. Setting ST to -2"
            )
            bestSet = self.bestAlleles(cg_sid)
            self.setPredictor(cg_sid, bestSet)
            return -2

    def bestST(self, cg_sid: str, st_list: list, type="profile"):
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
                    self.session.query(self.profiles[organism])
                    .filter(text(f"ST={st}"))
                    .first()
                )
            elif type == "novel":
                profiles.append(
                    self.session.query(self.novel[organism])
                    .filter(text(f"ST={st}"))
                    .first()
                )

        # Get values for each allele set that resolves an ST
        for prof in profiles:
            prof_keys = list(prof.keys())
            alleleconditions = list()
            alleledict = dict()
            allconditions = [f"Seq_types.CG_ID_sample=='{cg_sid}'"]

            for index, allele in enumerate(prof):
                if (
                    "ST" not in prof_keys[index]
                    and "clonal_complex" not in prof_keys[index]
                    and "species" not in prof_keys[index]
                ):
                    condition = f"Seq_types.loci=='{prof_keys[index]}' , Seq_types.allele=='{allele}'"
                    alleledict[prof_keys[index]] = ""
                    alleleconditions.append(f"and_({condition})")

            alleleconditions = f"or_({','.join(alleleconditions)})"
            allconditions.append(alleleconditions)
            allconditions = f"and_({','.join(allconditions)})"
            all_alleles = self.session.query(Seq_types).filter(eval(allconditions)).all()

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
                if not allele.loci in bestalleles[prof.ST].keys():
                    bestalleles[prof.ST][allele.loci] = dict()
                if not "contig_name" in bestalleles[prof.ST][allele.loci].keys():
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

    def bestAlleles(self, cg_sid: str):
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

    def get_unique_alleles(self, cg_sid: str, organism: str, threshold=True):
        """Returns a dict containing all unique alleles at every loci, and allele difference from expected"""
        tid = float(self.config["threshold"]["mlst_id"])
        tspan = (self.config["threshold"]["mlst_span"]) / 100.0
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
