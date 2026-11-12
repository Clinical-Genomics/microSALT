"""Compares existing organism references with available and updates as needed
By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
from typing import cast
import gzip
import json
import os
import re
import shutil
import shlex
import subprocess
import tempfile
import urllib.request
from urllib.parse import urlsplit
import xml.etree.ElementTree as ET
from pathlib import Path

from Bio import Entrez

from microSALT.config import (
    Folders,
    Threshold,
    PubMLSTCredentials,
    PasteurCredentials,
    Singularity,
    Containers,
)
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.exc.exceptions import MLSTUnavailableError, ReferenceSetupError
from microSALT.utils.pubmlst.client import BaseClient, get_client
from microSALT.utils.pubmlst.exceptions import InvalidURLError, PubMLSTError
from microSALT.utils.pubmlst.helpers import get_service_by_url


class Referencer:
    def __init__(
        self,
        log,
        folders: Folders,
        threshold: Threshold,
        pubmlst: PubMLSTCredentials,
        pasteur: PasteurCredentials,
        singularity: Singularity = None,
        containers: Containers = None,
        sampleinfo={},
        force=False,
    ):
        self.folders = folders
        self.threshold = threshold
        self.pubmlst = pubmlst
        self.pasteur = pasteur
        self.singularity = singularity or Singularity()
        self.containers = containers or Containers()
        self.logger = log
        self.db_access = DB_Manipulator(log=log, folders=folders, threshold=threshold)
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
        self.client: BaseClient = get_client(
            service, database, self.folders, self.pubmlst, self.pasteur
        )

    def _singularity_exec(self, tool: str, command: str) -> str:
        """Return command wrapped with singularity exec for the given tool container."""
        sif = getattr(self.containers, tool)
        binary = self.singularity.binary
        bind_list = list(self.singularity.bind_paths)
        bind = f"--bind {','.join(bind_list)}" if bind_list else ""
        return f"{binary} exec {bind} {sif} {command}"

    def identify_new(self) -> list[str]:
        """Download pubMLST & NCBI resources for organisms not already present.

        Returns a list of internal organism names (e.g. "staphylococcus_aureus")
        for which profile files were successfully downloaded. These names can be
        passed directly to create_new_profile_tables() to persist the profiles in
        the database without performing a full reference update.
        """
        neworgs: list[str] = []
        newrefs: list[str] = []
        downloaded: list[str] = []
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
                    f"{entry.get('reference')}.fasta" not in os.listdir(self.folders.genomes)
                    and entry.get("reference") not in newrefs
                ):
                    newrefs.append(cast(str, entry.get("reference")))
            for org in neworgs:
                truename = self.add_pubmlst(org)
                if truename is not None:
                    downloaded.append(truename)
            for org in newrefs:
                self.download_ncbi(org)
        except Exception:
            self.logger.error(
                "Unable to retrieve reference! Analysis using said reference will fail!"
            )
        return downloaded

    def create_new_profile_tables(self, new_organisms: list[str]) -> None:
        """Create database profile tables for organisms downloaded by identify_new().

        This is a lightweight alternative to update_refs() for the case where an
        organism is completely new to the system. Because the tables do not exist
        yet there is no risk of a concurrent reader seeing a partially-updated
        table, so no lock is acquired.
        """
        for organism in new_organisms:
            self.db_access.create_profile_table(organism)

    def prepare_references(self, samples, qc_only=False, skip_update=False) -> list[str]:
        """Check every required resource, collecting failures rather than hiding them."""
        errors = []
        references = {}
        organisms = {}
        needs_expec = False
        for sample in samples:
            name = str(sample.get("CG_ID_sample", "unknown sample"))
            reference = sample.get("reference")
            if (
                isinstance(reference, str)
                and reference
                and reference.lower() not in {"none", "null"}
            ):
                references.setdefault(reference, []).append(name)
            organism = sample.get("organism")
            if not qc_only and isinstance(organism, str) and organism:
                organisms.setdefault(organism, []).append(name)

        for reference, names in references.items():
            try:
                if not re.fullmatch(r"[\w.-]+", reference) or reference in {".", ".."}:
                    raise ReferenceSetupError("invalid reference identifier")
                path = Path(self.folders.genomes) / f"{reference}.fasta"
                if not path.exists():
                    if skip_update:
                        raise ReferenceSetupError(
                            f"missing reference {path}; downloads are disabled"
                        )
                    self.download_ncbi(reference)
                self.validate_fasta(path)
                missing = self._missing_genome_indexes(path)
                if missing:
                    if skip_update:
                        raise ReferenceSetupError(
                            f"missing or stale genome indexes: {', '.join(missing)}"
                        )
                    self._index_genome(path)
            except Exception as error:
                errors.append(f"Samples {', '.join(names)}, reference {reference}: {error}")

        for organism, names in organisms.items():
            try:
                ref = self.organism2reference(organism)
                needs_expec = needs_expec or ref == "escherichia_coli"
                if ref not in self.db_access.profiles:
                    if skip_update:
                        if ref is None:
                            self.logger.warning(
                                f"Samples {', '.join(names)}, organism {organism}: "
                                "no cached MLST scheme and downloads are disabled; skipping MLST"
                            )
                            continue
                        raise ReferenceSetupError("missing MLST profile; downloads are disabled")
                    ref = self.add_pubmlst(organism)
                    if not ref:
                        raise ReferenceSetupError("could not retrieve MLST resources")
                    self.create_new_profile_tables([ref])
                needs_expec = needs_expec or ref == "escherichia_coli"
                loci = Path(self.folders.references) / ref
                files = sorted(loci.glob("*.tfa"))
                if not files:
                    raise ReferenceSetupError(f"no MLST loci in {loci}")
                if not (Path(self.folders.profiles) / ref).is_file():
                    raise ReferenceSetupError(f"missing MLST profile file for {ref}")
                for path in files:
                    self.validate_fasta(path)
                self._prepare_blast_indexes(loci, ".tfa", skip_update)
            except MLSTUnavailableError as error:
                if ref is not None:
                    errors.append(
                        f"Samples {', '.join(names)}, organism {organism}: incomplete cached MLST resources; {error}"
                    )
                else:
                    self.logger.warning(
                        f"Samples {', '.join(names)}, organism {organism}: {error}; skipping MLST"
                    )
            except Exception as error:
                self.db_access.session.rollback()
                errors.append(f"Samples {', '.join(names)}, organism {organism}: {error}")

        if not qc_only:
            targets = [(Path(self.folders.resistances), ".fsa")]
            if needs_expec:
                targets.append((Path(self.folders.expec).parent, Path(self.folders.expec).suffix))
            for directory, suffix in targets:
                try:
                    files = sorted(directory.glob(f"*{suffix}"))
                    if not files:
                        raise ReferenceSetupError(f"no reference sequences in {directory}")
                    for path in files:
                        self.validate_fasta(path)
                    self._prepare_blast_indexes(directory, suffix, skip_update)
                except Exception as error:
                    errors.append(f"Typing resources {directory}: {error}")
        return errors

    def _prepare_blast_indexes(self, directory: Path, suffix: str, skip_update: bool):
        if not skip_update:
            self.index_db(str(directory), suffix)
        for path in directory.glob(f"*{suffix}"):
            for index in self._missing_blast_indexes(path):
                raise ReferenceSetupError(f"missing or stale BLAST index {index}")

    @staticmethod
    def _missing_blast_indexes(path: Path) -> list[Path]:
        indexes = [Path(f"{path.with_suffix('')}{suffix}") for suffix in (".nhr", ".nin", ".nsq")]
        return [
            index
            for index in indexes
            if not index.is_file()
            or not index.stat().st_size
            or index.stat().st_mtime < path.stat().st_mtime
        ]

    def update_refs(self):
        """Updates all references. Order is important, since no object is updated twice"""
        self.db_access.acquire_ref_lock()
        try:
            # Updates
            self.fetch_pubmlst(self.force)
            self.fetch_external()
            self.fetch_resistances(self.force)

            # Reindexes
            self.index_db(os.path.dirname(self.folders.expec), ".fsa")
        finally:
            self.db_access.release_ref_lock()

    def index_db(self, full_dir, suffix):
        """Rebuild missing/stale BLAST indexes and propagate indexing failures."""
        for path in sorted(Path(full_dir).glob(f"*{suffix}")):
            if not self._missing_blast_indexes(path):
                continue
            parse_ids = "" if suffix == ".fsa" else " -parse_seqids"
            command = self._singularity_exec(
                "blast",
                f"makeblastdb -in {shlex.quote(str(path))} -dbtype nucl{parse_ids} "
                f"-out {shlex.quote(path.stem)}",
            )
            self._run_index(command, full_dir)
            if self._missing_blast_indexes(path):
                raise ReferenceSetupError(
                    f"Indexing did not produce usable BLAST indexes for {path}"
                )

    def _parse_external_xml(
        self, url: str = "https://pubmlst.org/static/data/dbases.xml"
    ) -> ET.Element | None:
        """Fetch and parse the external XML, returning the root element."""
        try:
            query = urllib.request.urlopen(url).read()
            return ET.fromstring(query)
        except Exception as e:
            self.logger.warning(f"Unable to fetch or parse external XML: {e}")
            return None

    def _find_entry_for_organism(
        self, root: ET.Element, organism_name: str
    ) -> tuple[ET.Element | None, str | None]:
        """Find the XML entry for a given organism name."""
        organism_name = organism_name.lower().replace(" ", "_")
        for entry in root:
            species = entry.text.strip().lower().split(" ")[:2]
            organ: str = "_".join(species)
            if "escherichia_coli" in organ and "#1" in organ:
                if organism_name.startswith("escherichia"):
                    organ_no_suffix: str = organ[:-2]
                    return entry, organ_no_suffix
            if organ == organism_name:
                return entry, organ
        return None, None

    def _should_update_external(self, organ: str, entry: ET.Element) -> dict | bool:
        """Determine if the external data for an organism should be updated."""
        currver = self.db_access.read_version(f"profile_{organ}")
        st_link = entry.find("./mlst/database/profiles/url").text
        service = get_service_by_url(st_link)
        if service == "pasteur":
            database_url = entry.find(".//mlst/database/url").text.strip()
            if database_url is not None:
                database = database_url.split("/")[-2]
            database = f"pubmlst_{database}_seqdef"
            self.set_client(service, database=database)
        else:
            self.set_client(service)

        # Parse the database name and scheme ID
        try:
            parsed_data = self.client.parse_url(url=st_link)
        except InvalidURLError as e:
            self.logger.warning(f"Invalid URL: {st_link} - {e}")
            return False

        scheme_id = parsed_data.get("scheme_id")
        db = parsed_data.get("db")
        if not db or not scheme_id:
            self.logger.warning(
                f"Could not extract database name or scheme ID from MLST URL: {st_link}"
            )
            return False

        scheme_info = self.client.retrieve_scheme_info(db, scheme_id)
        last_updated = scheme_info.get("last_updated")
        if int(last_updated.replace("-", "")) <= int(currver.replace("-", "")) and not self.force:
            self.logger.info(
                f"Profile for {organ.replace('_', ' ').capitalize()} already at the latest version."
            )
            return False

        return {
            "organ": organ,
            "db": db,
            "scheme_id": scheme_id,
            "last_updated": last_updated,
            "scheme_info": scheme_info,
            "entry": entry,
            "service": service,
        }

    def _update_external_organism(self, update_info: dict) -> None:
        """Update the external data for a single organism."""
        organ = update_info["organ"]
        db = update_info["db"]
        scheme_id = update_info["scheme_id"]
        last_updated = update_info["last_updated"]
        scheme_info = update_info["scheme_info"]
        entry = update_info["entry"]
        service = update_info["service"]

        self.logger.info(
            f"{service} reference for {organ.replace('_', ' ').capitalize()} updated to {last_updated}"
        )

        # Step 1: Download the profiles CSV
        st_target = f"{self.folders.profiles}/{organ}"
        profiles_csv = self.client.download_profiles_csv(db, scheme_id)
        profiles_csv = profiles_csv.split("\n")
        trimmed_profiles = []
        trimmed_profiles.extend("\t".join(line.split("\t")[:8]) for line in profiles_csv)
        profiles_csv = "\n".join(trimmed_profiles)
        with open(st_target, "w") as profile_file:
            profile_file.write(profiles_csv)
        self.logger.info(f"Profiles CSV downloaded to {st_target}")

        # Step 2: Fetch scheme information to get loci
        loci_list = scheme_info.get("loci", [])

        # Step 3: Download loci FASTA files
        output = f"{self.folders.references}/{organ}"
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

        self.db_access.update_version(name=f"profile_{organ}", version=last_updated)
        self.db_access.reload_profiletable(organ)

    def fetch_external(self) -> None:
        """Fetch and update all external organisms as needed."""
        root = self._parse_external_xml()
        if root is None:
            return
        try:
            for entry in root:
                # Some species have extra names that are not expected, such as Klebsiella pneumoniae species complex, when we expect just Klebsiella pneumoniae
                species = entry.text.strip().lower().split(" ")[:2]
                organ = "_".join(species)
                # Escherichia coli and Acinetobacter baumannii have two entries, we use the one with the #1 suffix, but we have to remove it here
                if organ in ["escherichia_coli#1", "acinetobacter_baumannii#1"]:
                    organ = organ[:-2]
                if organ in self.organisms:
                    update_info = self._should_update_external(organ, entry)
                    if update_info:
                        self._update_external_organism(update_info)
        except PubMLSTError as e:
            self.logger.warning(f"Unable to update pubMLST external data: {e}")

    def fetch_external_for_organism(self, organism_name: str) -> None:
        """Fetch and update a specific organism from external sources."""
        root = self._parse_external_xml()
        if root is None:
            return
        entry, organ = self._find_entry_for_organism(root, organism_name)
        if entry and organ:
            update_info = self._should_update_external(organ, entry)
            if update_info:
                self._update_external_organism(update_info)
                self.logger.debug(f"Updated external data for {organism_name}")
        else:
            self.logger.warning(f"Organism '{organism_name}' not found in external sources.")

    def update_organism(self, external: bool, organism: str):
        """Fetches the latest information of an organism from the database"""
        if external:
            self.fetch_external_for_organism(organism)
        else:
            self.add_pubmlst(organism)

    def resync(self, type="", sample="", ignore=False):
        """Manipulates samples that have an internal ST that differs from pubMLST ST"""
        if type == "list":
            # Add single sample support later
            self.db_access.read_unresolved()
        elif type == "overwrite":
            if ignore:
                self.db_access.set_novel_ignored(sample=sample)
            else:
                self.db_access.set_novel_st(overwrite=True, sample=sample)
        else:
            self.db_access.set_novel_st(overwrite=False, sample=sample)

    def fetch_resistances(self, force=False):
        cwd = os.getcwd()
        url = "https://bitbucket.org/genomicepidemiology/resfinder_db.git"
        hiddensrc = f"{self.folders.resistances}/.resfinder_db"
        wipeIndex = False

        if not os.path.exists(hiddensrc) or len(os.listdir(hiddensrc)) == 0:
            self.logger.info("resFinder database not found. Caching..")
            if not os.path.exists(hiddensrc):
                os.makedirs(hiddensrc)
            cmd = f"git clone {url} --quiet"
            process = subprocess.Popen(
                cmd.split(),
                cwd=self.folders.resistances,
                stdout=subprocess.PIPE,
            )
            output, error = process.communicate()
            os.rename(
                f"{self.folders.resistances}/resfinder_db",
                hiddensrc,
            )
            wipeIndex = True
        else:
            if not wipeIndex:
                actual = os.listdir(self.folders.resistances)

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
                if "Already up-to-date." not in str(output):
                    self.logger.info("Remote resFinder database updated. Syncing...")
                    wipeIndex = True
                else:
                    self.logger.info("Cached resFinder database identical to remote.")

        # Actual update of resistance folder
        if wipeIndex:
            for file in os.listdir(hiddensrc):
                if os.path.isfile(f"{hiddensrc}/{file}"):
                    # Copy fresh
                    shutil.copy(
                        f"{hiddensrc}/{file}",
                        self.folders.resistances,
                    )

        # Double checks indexation is current.
        self.index_db(self.folders.resistances, ".fsa")

    def existing_organisms(self):
        """Returns list of all organisms currently added"""
        return self.organisms

    def organism2reference(self, normal_organism_name) -> str | None:
        """Finds which reference contains the same words as the organism
        and returns it in a format for database calls. Returns empty string if none found"""
        orgs = os.listdir(self.folders.references)
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
                f"Unable to find existing reference for {organism}, strain {normal_organism_name} has no reference match\nSource: {e}"
            )

    @staticmethod
    def validate_fasta(path: Path) -> None:
        """Reject error responses and malformed/empty nucleotide FASTA records."""
        records = 0
        length = 0
        with open(path) as stream:
            for number, line in enumerate(stream, 1):
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if not line[1:].strip() or (records and not length):
                        raise ReferenceSetupError(
                            f"{path}: empty FASTA header or sequence at line {number}"
                        )
                    records += 1
                    length = 0
                elif not records or not re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line):
                    raise ReferenceSetupError(
                        f"{path}: invalid nucleotide FASTA at line {number}: {line[:160]}"
                    )
                else:
                    length += len(line)
        if not records or not length:
            raise ReferenceSetupError(f"{path}: no non-empty FASTA sequences")

    @staticmethod
    def _missing_genome_indexes(path: Path) -> list[str]:
        missing = []
        for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"):
            index = Path(f"{path}{suffix}")
            if (
                not index.is_file()
                or not index.stat().st_size
                or index.stat().st_mtime < path.stat().st_mtime
            ):
                missing.append(str(index))
        return missing

    @staticmethod
    def _run_index(command: str, directory: str):
        try:
            subprocess.run(
                shlex.split(command), cwd=directory, check=True, capture_output=True, text=True
            )
        except subprocess.CalledProcessError as error:
            raise ReferenceSetupError(
                f"Index command failed (exit {error.returncode}): {command}: "
                f"{(error.stderr or error.stdout or '').strip()}"
            ) from error

    def _index_genome(self, path: Path):
        for tool, command in (("bwa", "bwa index"), ("samtools", "samtools faidx")):
            self._run_index(
                self._singularity_exec(tool, f"{command} {shlex.quote(str(path))}"),
                self.folders.genomes,
            )
        missing = self._missing_genome_indexes(path)
        if missing:
            raise ReferenceSetupError(
                f"Indexing did not produce usable files: {', '.join(missing)}"
            )

    def _assembly_fasta_url(self, reference: str) -> str:
        """Resolve the requested RefSeq/GenBank assembly, without substituting versions."""
        if not re.fullmatch(r"GC[AF]_\d+(?:\.\d+)?", reference):
            raise ReferenceSetupError(f"Invalid assembly accession {reference!r}")
        with Entrez.esearch(
            db="assembly", term=f"{reference}[Assembly Accession]", retmax=2, retmode="json"
        ) as handle:
            search = json.load(handle).get("esearchresult", {})
        identifiers = search.get("idlist", [])
        if int(search.get("count", 0)) != 1 or len(identifiers) != 1:
            raise ReferenceSetupError(
                f"Assembly accession {reference} did not resolve to exactly one NCBI Assembly record; "
                "check the accession and specify its version"
            )
        identifier = str(identifiers[0])
        with Entrez.esummary(db="assembly", id=identifier, retmode="json") as handle:
            summary = json.load(handle).get("result", {}).get(identifier, {})
        source = "refseq" if reference.startswith("GCF_") else "genbank"
        accession = summary.get("synonym", {}).get(source) or summary.get("assemblyaccession", "")
        matches = (
            accession == reference if "." in reference else accession.split(".")[0] == reference
        )
        if not matches:
            raise ReferenceSetupError(
                f"NCBI Assembly returned {accession!r} instead of requested accession {reference}"
            )
        location = summary.get(f"ftppath_{source}", "")
        if not location:
            raise ReferenceSetupError(
                f"No {source} genomic FASTA download is available for assembly {reference}"
            )
        parsed = urlsplit(location.rstrip("/"))
        directory = parsed.path.rsplit("/", 1)[-1]
        if (
            parsed.scheme not in {"ftp", "https"}
            or parsed.netloc != "ftp.ncbi.nlm.nih.gov"
            or not parsed.path.startswith("/genomes/all/")
            or parsed.query
            or parsed.fragment
            or not directory.startswith(f"{accession}_")
        ):
            raise ReferenceSetupError(
                f"Unexpected NCBI Assembly download path for {reference}: {location}"
            )
        return f"https://ftp.ncbi.nlm.nih.gov{parsed.path}/{directory}_genomic.fna.gz"

    def download_ncbi(self, reference):
        """Fetch a nucleotide or assembly genome, validate and index before publishing."""
        reference = "".join(reference.split())
        if not re.fullmatch(r"[\w.-]+", reference) or reference in {".", ".."}:
            raise ReferenceSetupError(f"Invalid reference identifier {reference!r}")
        try:
            Entrez.email = "2@2.com"
            with tempfile.TemporaryDirectory(
                prefix=".genome-", dir=self.folders.genomes
            ) as directory:
                staged = Path(directory) / f"{reference}.fasta"
                if reference.startswith(("GCF_", "GCA_")):
                    url = self._assembly_fasta_url(reference)
                    with (
                        urllib.request.urlopen(url, timeout=60) as response,
                        gzip.GzipFile(fileobj=response) as genome,
                        staged.open("wb") as target,
                    ):
                        shutil.copyfileobj(genome, target)
                else:
                    with (
                        Entrez.efetch(
                            db="nucleotide", id=reference, rettype="fasta", retmode="text"
                        ) as record,
                        staged.open("w") as target,
                    ):
                        shutil.copyfileobj(record, target)
                self.validate_fasta(staged)
                self._index_genome(staged)
                output = Path(self.folders.genomes) / staged.name
                # Publish the FASTA last: a failed fetch/index never replaces a cached genome.
                for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"):
                    Path(f"{staged}{suffix}").replace(f"{output}{suffix}")
                staged.replace(output)
            self.logger.info(f"Downloaded reference {reference}")
        except Exception as error:
            raise ReferenceSetupError(f"Unable to prepare genome '{reference}': {error}") from error

    def add_pubmlst(self, organism: str):
        """Checks pubmlst for references of given organism and downloads them"""
        # Organism must be in binomial format and only resolve to one hit
        errorg = organism
        self.db_access.check_ref_lock()
        try:
            organism = organism.lower().replace(".", " ")
            if organism.replace(" ", "_") in self.organisms and not self.force:
                self.logger.info(f"Organism {organism} already stored in microSALT")
                return organism.replace(" ", "_")
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
                            if part not in subtype["description"].lower():
                                missingPart = True
                    if not missingPart:
                        # Seqdef always appear after isolates, so this is fine
                        seqdef_url = subtype["href"]
                        desc = subtype["description"]
                        counter += 1.0
                        self.logger.info(f"Located pubMLST hit {desc} for sample")
            if counter > 2.0:
                raise Exception(
                    f"Reference '{errorg}' resolved to {int(counter / 2)} organisms. Please be more stringent"
                )
            elif counter < 1.0:
                raise MLSTUnavailableError(
                    f"No pubMLST database is available for organism '{errorg}'"
                )
            else:
                truename = desc.lower().split(" ")
                truename = f"{truename[0]}_{truename[1]}"
                self.download_pubmlst(truename, seqdef_url, force=self.force)
                # Update in-memory organism list so subsequent identify_new calls
                # in the same process see the new organism immediately.
                self.refs = self.db_access.profiles
                self.logger.info(f"Downloaded profile files for {truename}")
                return truename
        except MLSTUnavailableError:
            raise
        except Exception as e:
            raise ReferenceSetupError(f"Unable to prepare MLST resources for {errorg}: {e}") from e

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
                raise ReferenceSetupError(
                    f"Could not extract database name from URL: {subtype_href}"
                )

            # First, check scheme 1
            scheme_query_1 = self.client.retrieve_scheme_info(db, 1)
            mlst = None
            if "MLST" in (scheme_query_1 or {}).get("description", ""):
                mlst = f"{subtype_href}/schemes/1"
            else:
                # If scheme 1 isn't MLST, list all schemes and find the one with 'description' == 'MLST'
                record_query = self.client.list_schemes(db)
                for scheme in (record_query or {}).get("schemes", []):
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
            raise ReferenceSetupError(
                f"Unable to look up MLST schemes at {subtype_href}: {e}"
            ) from e

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
                raise ReferenceSetupError(
                    f"Could not extract database name or scheme ID from MLST URL: {mlst_href}"
                )

            scheme_info = self.client.retrieve_scheme_info(db, scheme_id)
            if scheme_info is None:
                raise MLSTUnavailableError(f"MLST scheme is no longer available at {mlst_href}")
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
        except MLSTUnavailableError:
            raise
        except Exception as e:
            raise ReferenceSetupError(
                f"Unable to determine pubMLST version for {organism}: {e}"
            ) from e

    def download_pubmlst(self, organism, subtype_href, force=False):
        """Download an available MLST scheme; return its version date when supplied."""
        organism = organism.lower().replace(" ", "_")
        try:
            # Pull version
            extver = self.external_version(organism, subtype_href)
            if not extver:
                self.logger.warning(
                    f"No external MLST version date for {organism}; skipping version comparison"
                )
            currver = self.db_access.read_version(f"profile_{organism}")
            if (
                extver
                and int(extver.replace("-", "")) <= int(currver.replace("-", ""))
                and not force
            ):
                self.logger.info(
                    f"Profile for {organism.replace('_', ' ').capitalize()} already at the latest version."
                )
                return currver

            # Retrieve the MLST scheme URL
            mlst_href = self.get_mlst_scheme(subtype_href)
            if not mlst_href:
                raise MLSTUnavailableError(f"No MLST scheme is available at {subtype_href}")

            # Parse the database name and scheme ID
            parsed_data = self.client.parse_url(url=mlst_href)
            db = parsed_data.get("db")
            scheme_id = parsed_data.get("scheme_id")
            if not db or not scheme_id:
                raise ReferenceSetupError(
                    f"Could not extract database name or scheme ID from MLST URL: {mlst_href}"
                )

            # Step 1: Download the profiles CSV
            st_target = f"{self.folders.profiles}/{organism}"
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
            if scheme_info is None:
                raise MLSTUnavailableError(f"MLST scheme is no longer available at {mlst_href}")
            loci_list = scheme_info.get("loci", [])

            # Step 3: Download loci FASTA files
            output = f"{self.folders.references}/{organism}"
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
        except MLSTUnavailableError:
            raise
        except Exception as e:
            raise ReferenceSetupError(f"Failed to download data for {organism}: {e}") from e

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
            internal_ver = self.db_access.read_version(f"profile_{key}")
            try:
                external_ver = self.external_version(key, val)

                if external_ver is None and not force:
                    self.logger.warning(
                        f"No external MLST version date for {key}; keeping cached resources"
                    )
                    continue
                if force or (internal_ver < external_ver):
                    self.logger.info(
                        f"pubMLST reference for {key.replace('_', ' ').capitalize()} updated to {external_ver} from {internal_ver}"
                    )
                    self.download_pubmlst(key, val, force)
                    if external_ver is not None:
                        self.db_access.update_version(name=f"profile_{key}", version=external_ver)
                    self.db_access.refresh_profiletable(key)
            except MLSTUnavailableError as error:
                self.logger.warning(f"{error}; keeping cached MLST resources for {key}")
