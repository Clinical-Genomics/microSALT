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
import zipfile

from Bio import Entrez
from bs4 import BeautifulSoup
from microSALT.store.db_manipulator import DB_Manipulator


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

    def identify_new(self, cg_id="", project=False):
        """ Automatically downloads pubMLST & NCBI organisms not already downloaded """
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
                "Reference update function failed prematurely. Review immediately"
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
        sufx_files = glob.glob(
            "{}/*{}".format(full_dir, suffix)
        )  # List of source files
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
                    if (
                        os.stat(file).st_mtime
                        < os.stat("{}/{}".format(full_dir, elem)).st_mtime
                    ):
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
                        bash_cmd = "makeblastdb -in {}/{} -dbtype nucl -parse_seqids -out {}".format(
                            full_dir, os.path.basename(file), os.path.basename(base)
                        )
                    proc = subprocess.Popen(
                        bash_cmd.split(), cwd=full_dir, stdout=subprocess.PIPE
                    )
                    output, error = proc.communicate()
                except Exception as e:
                    self.logger.error(
                        "Unable to index requested target {} in {}".format(
                            file, full_dir
                        )
                    )
        if reindexation:
            self.logger.info("Re-indexed contents of {}".format(full_dir))

    def fetch_external(self, force=False):
        """ Updates reference for data that IS ONLY LINKED to pubMLST """
        prefix = "https://pubmlst.org"
        query = urllib.request.urlopen("{}/data/".format(prefix))
        soup = BeautifulSoup(query, "html.parser")
        tr_sub = soup.find_all("tr", class_="td1")
        tr_sub = tr_sub + soup.find_all("tr", class_="td2")

        # Only search every other instance
        iterator = iter(tr_sub)
        unfound = True
        try:
            while unfound:
                entry = iterator.__next__()
                # Gather general info from first object
                sample = entry.get_text().split("\n")
                organ = sample[1].lower().replace(" ", "_")
                # In order to get ecoli #1
                if "escherichia_coli" in organ and "#1" in organ:
                    organ = organ[:-2]
                currver = self.db_access.get_version("profile_{}".format(organ))
                profile_no = re.search(r"\d+", sample[2]).group(0)
                if (
                    organ in self.organisms
                    and organ.replace("_", " ") not in self.updated
                    and (
                        int(profile_no.replace("-", "")) > int(currver.replace("-", ""))
                        or force
                    )
                ):
                    # Download definition files
                    st_link = prefix + entry.find_all("a")[1]["href"]
                    output = "{}/{}".format(self.config["folders"]["profiles"], organ)
                    urllib.request.urlretrieve(st_link, output)
                    # Update database
                    self.db_access.upd_rec(
                        {"name": "profile_{}".format(organ)},
                        "Versions",
                        {"version": profile_no},
                    )
                    self.db_access.reload_profiletable(organ)
                    # Gather loci from second object
                    entry = iterator.__next__()
                    # Clear existing directory and download allele files
                    out = "{}/{}".format(self.config["folders"]["references"], organ)
                    shutil.rmtree(out)
                    os.makedirs(out)
                    for loci in entry.find_all("a"):
                        loci = loci["href"]
                        lociname = os.path.basename(os.path.normpath(loci))
                        input = prefix + loci
                        urllib.request.urlretrieve(input, "{}/{}".format(out, lociname))
                    # Create new indexes
                    self.index_db(out, ".tfa")
                else:
                    iterator.__next__()
        except StopIteration:
            pass

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
                        self.logger.info(
                            "resFinder database files corrupted. Syncing..."
                        )
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
        """ Returns list of all organisms currently added """
        return self.organisms

    def organism2reference(self, normal_organism_name):
        """Finds which reference contains the same words as the organism
       and returns it in a format for database calls."""
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
            self.logger.warn(
                "Unable to find existing reference for {}, strain {} has no reference match\nSource: {}".format(
                    organism, normal_organism_name, e
                )
            )

    def download_ncbi(self, reference):
        """ Checks available references, downloads from NCBI if not present """
        try:
            DEVNULL = open(os.devnull, "wb")
            Entrez.email = "2@2.com"
            record = Entrez.efetch(
                db="nucleotide", id=reference, rettype="fasta", retmod="text"
            )
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
            self.logger.warning(
                "Unable to download genome '{}' from NCBI".format(reference)
            )

    def add_pubmlst(self, organism):
        """ Checks pubmlst for references of given organism and downloads them """
        # Organism must be in binomial format and only resolve to one hit
        errorg = organism
        try:
            organism = organism.lower().replace(".", " ")
            if organism.replace(" ", "_") in self.organisms and not self.force:
                self.logger.info(
                    "Organism {} already stored in microSALT".format(organism)
                )
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
                        self.logger.info(
                            "Located pubMLST hit {} for sample".format(desc)
                        )
            if counter > 2.0:
                raise Exception(
                    "Reference '{}' resolved to {} organisms. Please be more stringent".format(
                        errorg, int(counter / 2)
                    )
                )
            elif counter < 1.0:
                # add external
                raise Exception(
                    "Unable to find requested organism '{}' in pubMLST database".format(
                        errorg
                    )
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
        """ Returns a json object containing all organisms available via pubmlst.org """
        # Example request URI: http://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1/profiles_csv
        seqdef_url = dict()
        databases = "http://rest.pubmlst.org/db"
        db_req = urllib.request.Request(databases)
        with urllib.request.urlopen(db_req) as response:
            db_query = json.loads(response.read().decode("utf-8"))
        return db_query

    def download_pubmlst(self, organism, subtype_href, force=False):
        """ Downloads ST and loci for a given organism stored on pubMLST if it is more recent. Returns update date """
        organism = organism.lower().replace(" ", "_")

        # Pull version
        ver_req = urllib.request.Request("{}/schemes/1/profiles".format(subtype_href))
        with urllib.request.urlopen(ver_req) as response:
            ver_query = json.loads(response.read().decode("utf-8"))
        currver = self.db_access.get_version("profile_{}".format(organism))
        if (
            int(ver_query["last_updated"].replace("-", ""))
            <= int(currver.replace("-", ""))
            and not force
        ):
            # self.logger.info("Profile for {} already at latest version".format(organism.replace('_' ,' ').capitalize()))
            return currver

        # Pull ST file
        st_target = "{}/{}".format(self.config["folders"]["profiles"], organism)
        input = "{}/schemes/1/profiles_csv".format(subtype_href)
        urllib.request.urlretrieve(input, st_target)
        # Pull locus files
        loci_input = "{}/schemes/1".format(subtype_href)
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
        # Create new indexes
        self.index_db(output, ".tfa")

    def external_version(self, organism, subtype_href):
        ver_req = urllib.request.Request("{}/schemes/1/profiles".format(subtype_href))
        with urllib.request.urlopen(ver_req) as response:
            ver_query = json.loads(response.read().decode("utf-8"))
        return ver_query["last_updated"]

    def fetch_pubmlst(self, force=False):
        """ Updates reference for data that is stored on pubMLST """
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
            if internal_ver < external_ver:
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
