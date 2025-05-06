"""Compares existing organism references with available and updates as needed
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import glob
import os
import re
import shutil
import subprocess
import urllib.request

from microSALT.utils.pubmlst.client import BaseClient, PubMLSTClient, get_client
from microSALT.utils.pubmlst.authentication import ClientAuthentication

from Bio import Entrez
import xml.etree.ElementTree as ET
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.pubmlst.exceptions import InvalidURLError, PubMLSTError
from microSALT.utils.pubmlst.helpers import get_service_by_url


class Referencer:
    def __init__(self, config, log, sampleinfo={}, force=False):
        self.config = config
        self.logger = log
        self.db_access = DB_Manipulator(config, log)
        self.updated = list()
        # Fetch names of existing refs
        self.refs = self.db_access.profiles
        organisms = self.refs.keys()
        self.organisms = [*organisms]
        self.force = force

        self.sampleinfo = sampleinfo
        self.sample = None
        if isinstance(self.sampleinfo, list) and len(self.sampleinfo) > 1:
            self.name = self.sampleinfo[0].get("CG_ID_project")
            self.sample = self.sampleinfo[0]
            for entry in self.sampleinfo:
                if entry.get("CG_ID_sample") == self.name:
                    raise Exception(
                        "Mixed projects in samples_info file. Do not know how to proceed"
                    )
        else:
            if isinstance(self.sampleinfo, list):
                self.sampleinfo = self.sampleinfo[0]
            self.name = self.sampleinfo.get("CG_ID_sample")
            self.sample = self.sampleinfo
        self.client = None

    def set_client(self, service: str, database: str = None):
        """Set the client for PubMLST API interactions."""
        self.client: BaseClient = get_client(service, database)

    def identify_new(self, cg_id="", project=False):
        """Automatically downloads pubMLST & NCBI organisms not already downloaded"""
        neworgs = list()
        newrefs = list()
        try:
            if not isinstance(self.sampleinfo, list):
                samples = [self.sampleinfo]
            else:
                samples = self.sampleinfo

            for entry in samples:
                org = entry.get("organism")
                ref = self.organism2reference(org)
                if ref not in self.organisms and org not in neworgs:
                    neworgs.append(org)
                if (
                    not "{}.fasta".format(entry.get("reference"))
                    in os.listdir(self.config["folders"]["genomes"])
                    and not entry.get("reference") in newrefs
                ):
                    newrefs.append(entry.get("reference"))
            for org in neworgs:
                self.add_pubmlst(org)
            for org in newrefs:
                self.download_ncbi(org)
        except Exception as e:
            self.logger.error(
                "Unable to retrieve reference! Analysis using said reference will fail!"
            )

    def update_refs(self):
        """Updates all references. Order is important, since no object is updated twice"""
        # Updates
        self.fetch_pubmlst(self.force)
        self.fetch_external(self.force)
        self.fetch_resistances(self.force)

        # Reindexes
        self.index_db(os.path.dirname(self.config["folders"]["expec"]), ".fsa")

    def index_db(self, full_dir, suffix):
        """Check for indexation, makeblastdb job if not enough of them."""
        reindexation = False
        files = os.listdir(full_dir)
        sufx_files = glob.glob("{}/*{}".format(full_dir, suffix))  # List of source files
        for file in sufx_files:
            subsuf = "\{}$".format(suffix)
            base = re.sub(subsuf, "", file)

            bases = 0
            newer = 0
            for elem in files:
                # Number of files with same base name (7)
                if os.path.basename(base) == elem[: elem.rfind(".")]:
                    bases = bases + 1
                    # Number of index files fresher than source (6)
                    if os.stat(file).st_mtime < os.stat("{}/{}".format(full_dir, elem)).st_mtime:
                        newer = newer + 1
            # 7 for parse_seqids, 4 for not.
            if not (bases == 7 or newer == 6) and not (bases == 4 and newer == 3):
                reindexation = True
                try:
                    # Resistence files
                    if ".fsa" in suffix:
                        bash_cmd = "makeblastdb -in {}/{} -dbtype nucl -out {}".format(
                            full_dir, os.path.basename(file), os.path.basename(base)
                        )
                    # MLST locis
                    else:
                        bash_cmd = (
                            "makeblastdb -in {}/{} -dbtype nucl -parse_seqids -out {}".format(
                                full_dir, os.path.basename(file), os.path.basename(base)
                            )
                        )
                    proc = subprocess.Popen(bash_cmd.split(), cwd=full_dir, stdout=subprocess.PIPE)
                    proc.communicate()
                except Exception as e:
                    self.logger.error(
                        "Unable to index requested target {} in {}".format(file, full_dir)
                    )
        if reindexation:
            self.logger.info("Re-indexed contents of {}".format(full_dir))

    def fetch_external(self, force=False):
        url = "https://pubmlst.org/static/data/dbases.xml"
        try:
            query = urllib.request.urlopen(url).read()
            root = ET.fromstring(query)
            try:
                for entry in root:
                    # Check organism
                    species = entry.text.strip()
                    organ = species.lower().replace(" ", "_")
                    if "escherichia_coli" in organ and "#1" in organ:
                        organ = organ[:-2]
                    if organ in self.organisms:
                        # Check for newer version
                        currver = self.db_access.get_version("profile_{}".format(organ))
                        st_link = entry.find("./mlst/database/profiles/url").text
                        service: str = get_service_by_url(st_link)
                        if service == "pasteur":
                            database: str = f"pubmlst_{organ.split('_')[0]}_seqdef"
                            self.set_client(service, database=database)
                        else:
                            self.set_client(service)

                        # Parse the database name and scheme ID
                        try:
                            parsed_data = self.client.parse_url(url=st_link)
                        except InvalidURLError as e:
                            self.logger.warning(f"Invalid URL: {st_link} - {e}")
                            continue

                        scheme_id = parsed_data.get("scheme_id")  # Extract scheme ID
                        db = parsed_data.get("db")  # Extract database name

                        if not db or not scheme_id:
                            self.logger.warning(
                                f"Could not extract database name or scheme ID from MLST URL: {st_link}"
                            )
                            return

                        scheme_info = self.client.retrieve_scheme_info(
                            db, scheme_id
                        )  # Retrieve scheme info
                        last_updated = scheme_info.get("last_updated")  # Extract last updated date
                        if (
                            int(last_updated.replace("-", "")) <= int(currver.replace("-", ""))
                            and not force
                        ):
                            self.logger.info(
                                f"Profile for {organ.replace('_', ' ').capitalize()} already at the latest version."
                            )
                            continue
                        self.logger.info(
                            f"pubMLST reference for {organ.replace('_', ' ').capitalize()} updated to {last_updated} from {currver}"
                        )

                        # Step 1: Download the profiles CSV
                        st_target = f"{self.config['folders']['profiles']}/{organ}"
                        profiles_csv = self.client.download_profiles_csv(db, scheme_id)

                        # Only write the first 8 columns, this avoids adding information such as "clonal_complex" and "species"
                        profiles_csv = profiles_csv.split("\n")
                        trimmed_profiles = []
                        for line in profiles_csv:
                            trimmed_profiles.append("\t".join(line.split("\t")[:8]))

                        profiles_csv = "\n".join(trimmed_profiles)

                        with open(st_target, "w") as profile_file:
                            profile_file.write(profiles_csv)

                        self.logger.info(f"Profiles CSV downloaded to {st_target}")

                        # Step 2: Fetch scheme information to get loci

                        loci_list = scheme_info.get("loci", [])

                        # Step 3: Download loci FASTA files
                        output = f"{self.config['folders']['references']}/{organ}"
                        if os.path.isdir(output):
                            shutil.rmtree(output)
                        os.makedirs(output)

                        for locus_uri in loci_list:
                            locus_name = os.path.basename(os.path.normpath(locus_uri))
                            loci_fasta = self.client.download_locus(db, locus_name)
                            with open(f"{output}/{locus_name}.tfa", "w") as fasta_file:
                                fasta_file.write(loci_fasta)
                            self.logger.info(f"Locus FASTA downloaded: {locus_name}.tfa")

                        # Step 4: Create new indexes
                        self.index_db(output, ".tfa")

                        self.db_access.upd_rec(
                            {"name": "profile_{}".format(organ)},
                            "Versions",
                            {"version": last_updated},
                        )
                        self.db_access.reload_profiletable(organ)
            except PubMLSTError as e:
                self.logger.warning(f"Unable to update pubMLST external data: {e}")
        except Exception as e:
            self.logger.warning("Unable to update pubMLST external data: {}".format(e))

    def resync(self, type="", sample="", ignore=False):
        """Manipulates samples that have an internal ST that differs from pubMLST ST"""
        if type == "list":
            # Add single sample support later
            self.db_access.list_unresolved()
        elif type == "overwrite":
            if ignore:
                self.db_access.rm_novel(sample=sample)
            else:
                self.db_access.sync_novel(overwrite=True, sample=sample)
        else:
            self.db_access.sync_novel(overwrite=False, sample=sample)

    def fetch_resistances(self, force=False):
        cwd = os.getcwd()
        url = "https://bitbucket.org/genomicepidemiology/resfinder_db.git"
        hiddensrc = "{}/.resfinder_db".format(self.config["folders"]["resistances"])
        wipeIndex = False

        if not os.path.exists(hiddensrc) or len(os.listdir(hiddensrc)) == 0:
            self.logger.info("resFinder database not found. Caching..")
            if not os.path.exists(hiddensrc):
                os.makedirs(hiddensrc)
            cmd = "git clone {} --quiet".format(url)
            process = subprocess.Popen(
                cmd.split(),
                cwd=self.config["folders"]["resistances"],
                stdout=subprocess.PIPE,
            )
            output, error = process.communicate()
            os.rename(
                "{}/resfinder_db".format(self.config["folders"]["resistances"]),
                hiddensrc,
            )
            wipeIndex = True
        else:
            if not wipeIndex:
                actual = os.listdir(self.config["folders"]["resistances"])

                for file in os.listdir(hiddensrc):
                    if file not in actual and (".fsa" in file):
                        self.logger.info("resFinder database files corrupted. Syncing...")
                        wipeIndex = True
                        break

                cmd = "git pull origin master"
                process = subprocess.Popen(
                    cmd.split(),
                    cwd=hiddensrc,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )
                output, error = process.communicate()
                if not "Already up-to-date." in str(output):
                    self.logger.info("Remote resFinder database updated. Syncing...")
                    wipeIndex = True
                else:
                    self.logger.info("Cached resFinder database identical to remote.")

        # Actual update of resistance folder
        if wipeIndex:
            for file in os.listdir(hiddensrc):
                if os.path.isfile("{}/{}".format(hiddensrc, file)):
                    # Copy fresh
                    shutil.copy(
                        "{}/{}".format(hiddensrc, file),
                        self.config["folders"]["resistances"],
                    )

        # Double checks indexation is current.
        self.index_db(self.config["folders"]["resistances"], ".fsa")

    def existing_organisms(self):
        """Returns list of all organisms currently added"""
        return self.organisms

    def organism2reference(self, normal_organism_name):
        """Finds which reference contains the same words as the organism
        and returns it in a format for database calls. Returns empty string if none found"""
        orgs = os.listdir(self.config["folders"]["references"])
        organism = re.split(r"\W+", normal_organism_name.lower())
        try:
            for target in orgs:
                hit = 0
                for piece in organism:
                    if len(piece) == 1:
                        if target.startswith(piece):
                            hit += 1
                    elif piece in target or piece == "pneumonsiae" and "pneumoniae" in target:
                        hit += 1
                    else:
                        break
                if hit == len(organism):
                    return target
        except Exception as e:
            self.logger.warning(
                "Unable to find existing reference for {}, strain {} has no reference match\nSource: {}".format(
                    organism, normal_organism_name, e
                )
            )

    def download_ncbi(self, reference):
        """Checks available references, downloads from NCBI if not present"""
        try:
            DEVNULL = open(os.devnull, "wb")
            Entrez.email = "2@2.com"
            record = Entrez.efetch(db="nucleotide", id=reference, rettype="fasta", retmod="text")
            sequence = record.read()
            output = "{}/{}.fasta".format(self.config["folders"]["genomes"], reference)
            with open(output, "w") as f:
                f.write(sequence)
            bwaindex = "bwa index {}".format(output)
            proc = subprocess.Popen(
                bwaindex.split(),
                cwd=self.config["folders"]["genomes"],
                stdout=DEVNULL,
                stderr=DEVNULL,
            )
            out, err = proc.communicate()
            samindex = "samtools faidx {}".format(output)
            proc = subprocess.Popen(
                samindex.split(),
                cwd=self.config["folders"]["genomes"],
                stdout=DEVNULL,
                stderr=DEVNULL,
            )
            out, err = proc.communicate()
            self.logger.info("Downloaded reference {}".format(reference))
        except Exception as e:
            self.logger.warning("Unable to download genome '{}' from NCBI".format(reference))

    def add_pubmlst(self, organism):
        """Checks pubmlst for references of given organism and downloads them"""
        # Organism must be in binomial format and only resolve to one hit
        errorg = organism
        try:
            organism = organism.lower().replace(".", " ")
            if organism.replace(" ", "_") in self.organisms and not self.force:
                self.logger.info("Organism {} already stored in microSALT".format(organism))
                return
            db_query = self.query_pubmlst()

            # Doublecheck organism name is correct and unique
            orgparts = organism.split(" ")
            counter = 0.0
            for item in db_query:
                for subtype in item["databases"]:
                    missingPart = False
                    for part in orgparts:
                        if len(part) == 1:
                            if not subtype["description"].lower().startswith(part):
                                missingPart = True
                        else:
                            if not part in subtype["description"].lower():
                                missingPart = True
                    if not missingPart:
                        # Seqdef always appear after isolates, so this is fine
                        seqdef_url = subtype["href"]
                        desc = subtype["description"]
                        counter += 1.0
                        self.logger.info("Located pubMLST hit {} for sample".format(desc))
            if counter > 2.0:
                raise Exception(
                    "Reference '{}' resolved to {} organisms. Please be more stringent".format(
                        errorg, int(counter / 2)
                    )
                )
            elif counter < 1.0:
                # add external
                raise Exception(
                    "Unable to find requested organism '{}' in pubMLST database".format(errorg)
                )
            else:
                truename = desc.lower().split(" ")
                truename = "{}_{}".format(truename[0], truename[1])
                self.download_pubmlst(truename, seqdef_url)
                # Update organism list
                self.refs = self.db_access.profiles
                self.logger.info("Created table profile_{}".format(truename))
        except Exception as e:
            self.logger.warning(e.args[0])

    def query_pubmlst(self):
        """Returns a json object containing all organisms available via pubmlst.org"""
        self.set_client("pubmlst")
        db_query = self.client.query_databases()
        return db_query

    def get_mlst_scheme(self, subtype_href):
        """Returns the path for the MLST data scheme at pubMLST"""
        try:
            parsed_data = self.client.parse_url(url=subtype_href)
            db = parsed_data.get("db")
            if not db:
                self.logger.warning(f"Could not extract database name from URL: {subtype_href}")
                return None

            # First, check scheme 1
            scheme_query_1 = self.client.retrieve_scheme_info(db, 1)
            mlst = None
            if "MLST" in scheme_query_1.get("description", ""):
                mlst = f"{subtype_href}/schemes/1"
            else:
                # If scheme 1 isn't MLST, list all schemes and find the one with 'description' == 'MLST'
                record_query = self.client.list_schemes(db)
                for scheme in record_query.get("schemes", []):
                    if scheme.get("description") == "MLST":
                        mlst = scheme.get("scheme")
                        break

            if mlst:
                self.logger.debug(f"Found data at pubMLST: {mlst}")
                return mlst
            else:
                self.logger.warning(f"Could not find MLST data at {subtype_href}")
                return None
        except Exception as e:
            self.logger.warning(e)
            return None

    def external_version(self, organism, subtype_href):
        """Returns the version (date) of the data available on pubMLST"""
        try:
            mlst_href = self.get_mlst_scheme(subtype_href)
            if not mlst_href:
                self.logger.warning(f"MLST scheme not found for URL: {subtype_href}")
                return None
            parsed_data = self.client.parse_url(url=mlst_href)
            db = parsed_data.get("db")
            scheme_id = parsed_data.get("scheme_id")
            if not db or not scheme_id:
                self.logger.warning(
                    f"Could not extract database name or scheme ID from MLST URL: {mlst_href}"
                )
                return None

            scheme_info = self.client.retrieve_scheme_info(db, scheme_id)
            last_updated = scheme_info.get("last_updated")
            if last_updated:
                self.logger.debug(
                    f"Retrieved last_updated: {last_updated} for organism: {organism}"
                )
                return last_updated
            else:
                self.logger.warning(
                    f"No 'last_updated' field found for db: {db}, scheme_id: {scheme_id}"
                )
                return None
        except Exception as e:
            self.logger.warning(f"Could not determine pubMLST version for {organism}")
            self.logger.warning(e)
            return None

    def download_pubmlst(self, organism, subtype_href, force=False):
        """Downloads ST and loci for a given organism stored on pubMLST if it is more recent. Returns update date"""
        organism = organism.lower().replace(" ", "_")
        try:
            # Pull version
            extver = self.external_version(organism, subtype_href)
            currver = self.db_access.get_version(f"profile_{organism}")
            if int(extver.replace("-", "")) <= int(currver.replace("-", "")) and not force:
                self.logger.info(
                    f"Profile for {organism.replace('_', ' ').capitalize()} already at the latest version."
                )
                return currver

            # Retrieve the MLST scheme URL
            mlst_href = self.get_mlst_scheme(subtype_href)
            if not mlst_href:
                self.logger.warning(f"MLST scheme not found for URL: {subtype_href}")
                return None

            # Parse the database name and scheme ID
            parsed_data = self.client.parse_url(url=mlst_href)
            db = parsed_data.get("db")
            scheme_id = parsed_data.get("scheme_id")
            if not db or not scheme_id:
                self.logger.warning(
                    f"Could not extract database name or scheme ID from MLST URL: {mlst_href}"
                )
                return None

            # Step 1: Download the profiles CSV
            st_target = f"{self.config['folders']['profiles']}/{organism}"
            profiles_csv = self.client.download_profiles_csv(db, scheme_id)
            # Only write the first 8 columns, this avoids adding information such as "clonal_complex" and "species"
            profiles_csv = profiles_csv.split("\n")
            trimmed_profiles = []
            for line in profiles_csv:
                trimmed_profiles.append("\t".join(line.split("\t")[:8]))

            profiles_csv = "\n".join(trimmed_profiles)

            with open(st_target, "w") as profile_file:
                profile_file.write(profiles_csv)
            self.logger.info(f"Profiles CSV downloaded to {st_target}")

            # Step 2: Fetch scheme information to get loci
            scheme_info = self.client.retrieve_scheme_info(db, scheme_id)
            loci_list = scheme_info.get("loci", [])

            # Step 3: Download loci FASTA files
            output = f"{self.config['folders']['references']}/{organism}"
            if os.path.isdir(output):
                shutil.rmtree(output)
            os.makedirs(output)

            for locus_uri in loci_list:
                locus_name = os.path.basename(os.path.normpath(locus_uri))
                loci_fasta = self.client.download_locus(db, locus_name)
                with open(f"{output}/{locus_name}.tfa", "w") as fasta_file:
                    fasta_file.write(loci_fasta)
                self.logger.info(f"Locus FASTA downloaded: {locus_name}.tfa")

            # Step 4: Create new indexes
            self.index_db(output, ".tfa")

            return extver
        except Exception as e:
            self.logger.error(f"Failed to download data for {organism}: {e}")
            return None

    def fetch_pubmlst(self, force=False):
        """Updates reference for data that is stored on pubMLST"""
        seqdef_url = dict()
        db_query = self.query_pubmlst()

        # Fetch seqdef locations
        for item in db_query:
            for subtype in item["databases"]:
                for name in self.organisms:
                    if name.replace("_", " ") in subtype["description"].lower():
                        # Seqdef always appear after isolates, so this is fine
                        self.updated.append(name.replace("_", " "))
                        seqdef_url[name] = subtype["href"]

        for key, val in seqdef_url.items():
            internal_ver = self.db_access.get_version("profile_{}".format(key))
            external_ver = self.external_version(key, val)

            if (internal_ver < external_ver) or force:
                self.logger.info(
                    "pubMLST reference for {} updated to {} from {}".format(
                        key.replace("_", " ").capitalize(), external_ver, internal_ver
                    )
                )
                self.download_pubmlst(key, val, force)
                self.db_access.upd_rec(
                    {"name": "profile_{}".format(key)},
                    "Versions",
                    {"version": external_ver},
                )
                self.db_access.reload_profiletable(key)
