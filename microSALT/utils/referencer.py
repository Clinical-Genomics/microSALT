"""Compares existing organism references with available and updates as needed
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import glob
import json
import os
import re
import shutil
import subprocess
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

from Bio import Entrez

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.pubmlst.api import (
    check_database_metadata,
    download_locus,
    download_profiles_csv,
    fetch_schemes,
    query_databases,
)
from microSALT.utils.pubmlst.authentication import (
    get_new_session_token,
    load_session_token,
)

def resolve_path(path):
    """Resolve environment variables, user shortcuts, and convert to absolute path."""
    if path:
        path = os.path.expandvars(path)
        path = os.path.expanduser(path)
        path = os.path.abspath(path)
    return path

def resolve_config_paths(self, config):
    # Ensure all paths in 'folders' are resolved
    if 'folders' in config:
        for key, value in config['folders'].items():
            if isinstance(value, str) and '/' in value:
                config['folders'][key] = resolve_path(value)

    # Resolve pubmlst credentials_files_path if present
    if 'pubmlst' in config and 'credentials_files_path' in config['pubmlst']:
        config['pubmlst']['credentials_files_path'] = resolve_path(config['pubmlst']['credentials_files_path'])

    return config

def ensure_directories(self):
    """Ensure all required directories are created."""
    required_dirs = [
        self.config["folders"].get("results"),
        self.config["folders"].get("reports"),
        self.config["folders"].get("profiles"),
        self.config["folders"].get("references"),
        self.config["folders"].get("resistances"),
        self.config["folders"].get("genomes"),
    ]
    for dir_path in required_dirs:
        if dir_path:
            resolved_path = resolve_path(dir_path)
            os.makedirs(resolved_path, exist_ok=True)
            self.logger.info(f"Ensured directory exists: {resolved_path}")


class Referencer:
    def __init__(self, config, log, sampleinfo={}, force=False):
        self.config = self.resolve_config_paths(config)
        self.ensure_directories() 
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

        # Use a default database to load or fetch an initial token
        default_db = "pubmlst_test_seqdef"
        self.token, self.secret = load_session_token(default_db)
        if not self.token or not self.secret:
            self.token, self.secret = get_new_session_token(default_db)

    def resolve_config_paths(self, config):
        # Ensure all paths in 'folders' are resolved
        if 'folders' in config:
            for key, value in config['folders'].items():
                if isinstance(value, str) and '/' in value:
                    config['folders'][key] = resolve_path(value)

        # Resolve pubmlst credentials_files_path if present
        if 'pubmlst' in config and 'credentials_files_path' in config['pubmlst']:
            config['pubmlst']['credentials_files_path'] = resolve_path(config['pubmlst']['credentials_files_path'])

        return config

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
        full_dir = resolve_path(full_dir)
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
                    profiles_query = urllib.request.urlopen(st_link)
                    profile_no = profiles_query.readlines()[-1].decode("utf-8").split("\t")[0]
                    if organ.replace("_", " ") not in self.updated and (
                        int(profile_no.replace("-", "")) > int(currver.replace("-", "")) or force
                    ):
                        # Download MLST profiles
                        self.logger.info("Downloading new MLST profiles for " + species)
                        output = "{}/{}".format(self.config["folders"]["profiles"], organ)
                        urllib.request.urlretrieve(st_link, output)
                        # Clear existing directory and download allele files
                        out = "{}/{}".format(self.config["folders"]["references"], organ)
                        shutil.rmtree(out)
                        os.makedirs(out)
                        for locus in entry.findall("./mlst/database/loci/locus"):
                            locus_name = locus.text.strip()
                            locus_link = locus.find("./url").text
                            urllib.request.urlretrieve(
                                locus_link, "{}/{}.tfa".format(out, locus_name)
                            )
                        # Create new indexes
                        self.index_db(out, ".tfa")
                        # Update database
                        self.db_access.upd_rec(
                            {"name": "profile_{}".format(organ)},
                            "Versions",
                            {"version": profile_no},
                        )
                        self.db_access.reload_profiletable(organ)
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
            refs = 0
            for target in orgs:
                hit = 0
                for piece in organism:
                    if len(piece) == 1:
                        if target.startswith(piece):
                            hit += 1
                    else:
                        if piece in target:
                            hit += 1
                        # For when people misspell the strain in the orderform
                        elif piece == "pneumonsiae" and "pneumoniae" in target:
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
        errorg = organism
        try:
            organism = organism.lower().replace(".", " ")
            if organism.replace(" ", "_") in self.organisms and not self.force:
                self.logger.info("Organism {} already stored in microSALT".format(organism))
                return
            db_query = self.query_pubmlst()

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
                raise Exception(
                    "Unable to find requested organism '{}' in pubMLST database".format(errorg)
                )
            else:
                truename = desc.lower().split(" ")
                truename = "{}_{}".format(truename[0], truename[1])
                self.download_pubmlst(truename, seqdef_url)
                self.refs = self.db_access.profiles
                self.logger.info("Created table profile_{}".format(truename))
        except Exception as e:
            self.logger.warning(e.args[0])

    def query_pubmlst(self):
        """Returns a json object containing all organisms available via pubmlst.org"""
        databases = "http://rest.pubmlst.org/db"
        db_req = urllib.request.Request(databases)
        with urllib.request.urlopen(db_req) as response:
            db_query = json.loads(response.read().decode("utf-8"))
        return db_query

    def get_mlst_scheme(self, subtype_href):
        """Returns the path for the MLST data scheme at pubMLST"""
        try:
            mlst = False
            record_req_1 = urllib.request.Request("{}/schemes/1".format(subtype_href))
            with urllib.request.urlopen(record_req_1) as response:
                scheme_query_1 = json.loads(response.read().decode("utf-8"))
                if "MLST" in scheme_query_1["description"]:
                    mlst = "{}/schemes/1".format(subtype_href)
            if not mlst:
                record_req = urllib.request.Request("{}/schemes".format(subtype_href))
                with urllib.request.urlopen(record_req) as response:
                    record_query = json.loads(response.read().decode("utf-8"))
                    for scheme in record_query["schemes"]:
                        if scheme["description"] == "MLST":
                            mlst = scheme["scheme"]
            if mlst:
                self.logger.debug("Found data at pubMLST: {}".format(mlst))
                return mlst
            else:
                self.logger.warning("Could not find MLST data at {}".format(subtype_href))
        except Exception as e:
            self.logger.warning(e)

    def external_version(self, organism, subtype_href):
        """Returns the version (date) of the data available on pubMLST"""
        mlst_href = self.get_mlst_scheme(subtype_href)
        try:
            with urllib.request.urlopen(mlst_href) as response:
                ver_query = json.loads(response.read().decode("utf-8"))
            return ver_query["last_updated"]
        except Exception as e:
            self.logger.warning("Could not determine pubMLST version for {}".format(organism))
            self.logger.warning(e)

    def download_pubmlst(self, organism, subtype_href, force=False):
        """Downloads ST and loci for a given organism stored on pubMLST if it is more recent. Returns update date"""
        organism = organism.lower().replace(" ", "_")

        extver = self.external_version(organism, subtype_href)
        currver = self.db_access.get_version("profile_{}".format(organism))
        if int(extver.replace("-", "")) <= int(currver.replace("-", "")) and not force:
            return currver

        mlst_href = self.get_mlst_scheme(subtype_href)
        st_target = "{}/{}".format(self.config["folders"]["profiles"], organism)
        st_input = "{}/profiles_csv".format(mlst_href)
        urllib.request.urlretrieve(st_input, st_target)

        loci_input = mlst_href
        loci_req = urllib.request.Request(loci_input)
        with urllib.request.urlopen(loci_req) as response:
            loci_query = json.loads(response.read().decode("utf-8"))

        output = "{}/{}".format(self.config["folders"]["references"], organism)

        try:
            if os.path.isdir(output):
                shutil.rmtree(output)
        except FileNotFoundError as e:
            pass
        os.makedirs(output)

        for locipath in loci_query["loci"]:
            loci = os.path.basename(os.path.normpath(locipath))
            urllib.request.urlretrieve(
                "{}/alleles_fasta".format(locipath), "{}/{}.tfa".format(output, loci)
            )
        self.index_db(output, ".tfa")

    def fetch_pubmlst(self, force=False):
        """Fetches and updates PubMLST data."""
        try:
            self.logger.info("Querying available PubMLST databases...")
            databases = query_databases(self.token, self.secret)

            for db_entry in databases:
                db_name = db_entry["name"]
                db_desc = db_entry["description"]

                for sub_db in db_entry.get("databases", []):
                    sub_db_name = sub_db["name"]
                    sub_db_desc = sub_db["description"]

                    # Skip databases that are not sequence definitions or do not match known organisms
                    if "seqdef" not in sub_db_name.lower():
                        self.logger.debug(f"Skipping {sub_db_desc} (not a sequence definition database).")
                        continue

                    if sub_db_desc.replace(" ", "_").lower() not in self.organisms and not force:
                        self.logger.debug(f"Skipping {sub_db_desc}, not in known organisms.")
                        continue

                    # Load or fetch a session token for this specific sub-database
                    db_token, db_secret = load_session_token(sub_db_name)
                    if not db_token or not db_secret:
                        db_token, db_secret = get_new_session_token(sub_db_name)

                    self.logger.info(f"Fetching schemes for {sub_db_desc}...")
                    schemes = fetch_schemes(sub_db_name, db_token, db_secret)

                    for scheme in schemes.get("schemes", []):
                        if "scheme" not in scheme:
                            self.logger.warning(f"Scheme does not contain 'scheme' key: {scheme}")
                            continue

                        scheme_url = scheme["scheme"]
                        scheme_id = scheme_url.rstrip("/").split("/")[-1]

                        if not scheme_id.isdigit():
                            self.logger.error(f"Invalid scheme ID: {scheme_url}")
                            continue

                        if "MLST" in scheme["description"]:
                            self.logger.debug(f"Downloading profiles for {sub_db_desc}...")
                            # Use the CSV endpoint to avoid pagination issues
                            try:
                                profiles_csv = download_profiles_csv(sub_db_name, scheme_id, db_token, db_secret)
                                org_folder_name = sub_db_desc.replace(" ", "_").lower()
                                st_target = "{}/{}".format(self.config["folders"]["profiles"], org_folder_name)
                                with open(st_target, "w") as f:
                                    f.write(profiles_csv)

                                # Process loci
                                loci = scheme.get("loci", [])
                                if not loci:
                                    self.logger.warning(f"No loci found for scheme {scheme_id} in {sub_db_desc}.")
                                else:
                                    out = "{}/{}".format(self.config["folders"]["references"], org_folder_name)
                                    if os.path.isdir(out):
                                        shutil.rmtree(out)
                                    os.makedirs(out)
                                    for locus in loci:
                                        self.logger.info(f"Downloading locus {locus} for {sub_db_desc}...")
                                        locus_data = download_locus(sub_db_name, locus, db_token, db_secret)
                                        locus_file_path = os.path.join(out, f"{locus}.tfa")
                                        with open(locus_file_path, "wb") as locus_file:
                                            locus_file.write(locus_data)
                                        self.logger.info(f"Locus {locus} downloaded successfully.")
                                    self.index_db(out, ".tfa")

                                # Check and log metadata
                                metadata = check_database_metadata(sub_db_name, db_token, db_secret)
                                last_updated = metadata.get("last_updated", "Unknown")
                                if last_updated != "Unknown":
                                    self.db_access.upd_rec(
                                        {"name": f"profile_{org_folder_name}"},
                                        "Versions",
                                        {"version": last_updated},
                                    )
                                    self.db_access.reload_profiletable(org_folder_name)
                                    self.logger.info(f"Database {sub_db_desc} updated to {last_updated}.")
                                else:
                                    self.logger.debug(f"No new updates for {sub_db_desc}.")
                            except Exception as e:
                                self.logger.error(f"Error processing {sub_db_desc}: {e}")

            self.logger.info("PubMLST fetch and update process completed successfully.")
        except Exception as e:
            self.logger.error(f"Failed to fetch PubMLST data: {e}")
