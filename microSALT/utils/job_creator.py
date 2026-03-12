"""Creates sbatch jobs for MLST instances
By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import glob
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from importlib.metadata import entry_points
from pathlib import Path

import yaml

from microSALT import __version__
from microSALT.config import (
    Containers,
    Folders,
    PasteurCredentials,
    PubMLSTCredentials,
    Regex,
    Singularity,
    SlurmHeader,
    Threshold,
)
from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.utils.referencer import Referencer


class Job_Creator:
    def __init__(
        self,
        log,
        folders: Folders,
        slurm_header: SlurmHeader,
        regex: Regex,
        dry: bool,
        config_path: str,
        threshold: Threshold,
        pubmlst: PubMLSTCredentials,
        pasteur: PasteurCredentials,
        singularity: Singularity,
        containers: Containers,
        sampleinfo={},
        run_settings={},
    ):
        self.folders = folders
        self.slurm_header = slurm_header
        self.regex = regex
        self.dry = dry
        self.config_path = config_path
        self.threshold = threshold
        self.pubmlst = pubmlst
        self.pasteur = pasteur
        self.singularity = singularity
        self.containers = containers
        self.logger = log
        self.batchfile = "/tmp/batchfile.sbatch"

        self.filelist = list()
        if isinstance(run_settings.get("input"), list):
            self.filelist = run_settings.get("input")
            run_settings["input"] = "/tmp/"

        self.run_settings = run_settings
        self.indir = os.path.abspath(run_settings.get("input", "/tmp/"))
        self.trimmed = run_settings.get("trimmed", True)
        self.qc_only = run_settings.get("qc_only", False)
        self.pool = run_settings.get("pool", [])
        self.finishdir = run_settings.get("finishdir", "")

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

        # If timestamp is provided. Use it as analysis time. Else use current time
        if run_settings.get("timestamp") is not None:
            self.now = run_settings.get("timestamp")
            temp = run_settings.get("timestamp").replace("_", ".").split(".")
            self.dt = datetime(
                int(temp[0]),
                int(temp[1]),
                int(temp[2]),
                int(temp[3]),
                int(temp[4]),
                int(temp[5]),
            )
        else:
            self.dt = datetime.now()
            self.now = time.strftime(
                f"{self.dt.year}.{self.dt.month}.{self.dt.day}_{self.dt.hour}.{self.dt.minute}.{self.dt.second}"
            )

        if run_settings.get("finishdir") is None:
            self.finishdir = f"{folders.results}/{self.name}_{self.now}"
        self.db_pusher = DB_Manipulator(log=log, folders=folders, threshold=threshold)
        self.concat_files = dict()
        self.ref_resolver = Referencer(
            log=log,
            folders=folders,
            threshold=threshold,
            pubmlst=pubmlst,
            pasteur=pasteur,
            singularity=singularity,
            containers=containers,
        )

    def get_sbatch(self):
        """Returns sbatchfile, slightly superflous"""
        return self.batchfile

    def get_headerargs(self):
        headerline = f"-A {self.slurm_header.project} -p {self.slurm_header.type} -n {self.slurm_header.threads} -t {self.slurm_header.time} -J {self.slurm_header.job_prefix}_{self.name} --qos {self.slurm_header.qos} --output {self.finishdir}/slurm_{self.name}.log"
        return headerline

    def _singularity_exec(self, tool: str, command: str) -> str:
        """Return command wrapped with singularity exec for the given tool container."""
        sif = getattr(self.containers, tool)
        binary = self.singularity.binary
        bind_list = list(self.singularity.bind_paths)
        if self.finishdir and self.finishdir not in bind_list:
            bind_list.append(self.finishdir)
        bind = f"--bind {','.join(bind_list)}" if bind_list else ""
        return f"{binary} exec {bind} {sif} {command}"

    def verify_fastq(self):
        """Uses arg indir to return a dict of PE fastq tuples fulfilling naming convention"""
        verified_files = list()
        files = os.listdir(self.indir)
        if files == []:
            raise Exception(f"Directory {self.indir} lacks fastq files.")
        for file in files:
            file_match = re.match(self.regex.file_pattern, file)
            if file_match:
                # Check that symlinks resolve
                path = f"{self.indir}/{file}"
                if os.path.islink(path):
                    if not os.path.exists(os.readlink(path)):
                        raise Exception(
                            f"Some fastq files are unresolved symlinks in directory {self.indir}."
                        )

                # Make sure both mates exist
                if (
                    file_match[1] == "1"
                    or file_match[1] == "2"
                    or file_match[1] == "forward"
                    or file_match[1] == "reverse"
                ):
                    if file_match[1] == "forward" or file_match[1] == "reverse":
                        pairno = "forward"
                        if "forward" in file_match[1]:
                            pairno = "reverse"
                        pairname = file_match[0].replace(file_match[1], pairno)
                    else:
                        pairno = 2 - 1 % int(file_match[1])  # 1->2, 2->1
                        # Construct mate name
                        pairname = f"{file_match.string[: file_match.end(1) - 1]}{pairno}{file_match.string[file_match.end(1) : file_match.end()]}"
                    if pairname in files:
                        files.pop(files.index(pairname))
                        verified_files.append(file_match[0])
                        verified_files.append(pairname)
                else:
                    raise Exception(f"Some fastq files have no mate in directory {self.indir}.")
        if verified_files == []:
            raise Exception(
                f"No files in directory {self.indir} match file_pattern '{self.regex.file_pattern}'."
            )

        # Warn about file sizes
        for vfile in verified_files:
            try:
                bsize = os.stat(f"{self.indir}/{vfile}").st_size
                bsize = bsize >> 20
                if bsize > 1000:
                    self.logger.warning(f"Input fastq {vfile} exceeds 1000MB")
            except Exception:
                self.logger.warning(f"Unable to verify size of input file {self.indir}/{vfile}")

        # Warn about invalid fastq files
        for vfile in verified_files:
            f = gzip.open(f"{self.indir}/{vfile}", "r")
            lines = f.read().splitlines()
            if len(lines) < 2 or "+" not in str(lines[-2]):
                self.logger.warning(f"Input fastq {vfile} does not seem to end properly")
        return sorted(verified_files)

    @staticmethod
    def create_version_file(finishdir):
        """Creates a section in the batchfile with the version of microSALT"""
        version_file_path = f"{finishdir}/version.txt"
        with open(version_file_path, "w") as version_file:
            version_file.write(__version__)

    def create_assemblysection(self):
        assembly_dir = f"{self.finishdir}/assembly"
        contigs_file_raw = f"{assembly_dir}/{self.name}_contigs_raw.fasta"
        contigs_file = f"{assembly_dir}/{self.name}_contigs.fasta"
        contigs_trimmed_file = f"{assembly_dir}/{self.name}_trimmed_contigs.fasta"

        batchfile = open(self.batchfile, "a+")
        # memory is actually 128 per node regardless of cores.
        batchfile.write("# SKESA assembly\n")
        skesa_cmd = (
            f"skesa "
            f"--cores {self.slurm_header.threads} "
            f"--memory {8 * int(self.slurm_header.threads)} "
            f"--contigs_out {contigs_file_raw} "
            f"--reads {self.concat_files['f']},{self.concat_files['r']}"
        )
        batchfile.write(f"mkdir -p {assembly_dir} &{self._singularity_exec('skesa', skesa_cmd)}\n")

        # Convert sequence naming in Skesa output into Spades format in the contigs fasta file:
        # ----------------------------------------------
        # Skesa format:  >Contig_1_100.000
        # Spades format: >NODE_1_length_150_cov_100.000
        # ----------------------------------------------
        # We do the change by doing the following with awk:
        # 1. When the line is a header (starting with >), capture the contig number and coverage
        # 2. When the line is NOT a header (not starting with >), compute the length of the line
        #    and then print out:
        #    a. A new header line in Spades format with the length included
        #    b. Print out the sequence line.
        # Note: The match function requires GNU awk (gawk) to be able to capture groups in regexes.
        batchfile.write(
            "gawk "
            + r"'/^>/ { match($0, /Contig_([0-9]+)_([0-9\.]+)/, m) } "
            + '!/^>/ { seqlen=length($0); print ">NODE_" m[1] "_length_" seqlen "_cov_" m[2]; print $0; }\' '
            + f"{contigs_file_raw} > {contigs_file}\n"
        )

        # Keep only the 999(?) top contigs to avoid really low-quality contigs
        batchfile.write(f"sed -n '/NODE_1000_/q;p' {contigs_file} > {contigs_trimmed_file}\n")
        # batchfile.write("##Input cleanup\n")
        # batchfile.write("rm -r {}/trimmed\n".format(self.finishdir))
        batchfile.write("\n\n")
        batchfile.close()

    def blast_subset(self, name, search_string):
        # Create run
        file_list = glob.glob(search_string)
        batchfile = open(self.batchfile, "a+")
        batchfile.write(f"mkdir {self.finishdir}/blast_search/{name}\n")
        blast_format = '"7 stitle sstrand qaccver saccver pident evalue bitscore qstart qend sstart send length"'

        if len(file_list) > 1:
            for ref in file_list:
                if re.search(r"(\w+(?:\-\w+)*)\.\w+", os.path.basename(ref)) is None:
                    self.logger.error(
                        "File {} does not match typical format. Consider deleting and redownloading"
                    )
                else:
                    ref_nosuf = re.search(r"(\w+(?:\-\w+)*)\.\w+", os.path.basename(ref)).group(1)
                batchfile.write(
                    f"# BLAST {name} search for {self.sample.get('organism')}, {ref_nosuf}\n"
                )
                if name == "mlst":
                    batchfile.write(
                        self._singularity_exec(
                            "blast",
                            f"blastn -db {os.path.dirname(ref)}/{ref_nosuf}  -query {self.finishdir}/assembly/{self.name}_contigs.fasta -out {self.finishdir}/blast_search/{name}/loci_query_{ref_nosuf}.txt -task megablast -num_threads {self.slurm_header.threads} -outfmt {blast_format}",
                        )
                        + "\n"
                    )
                else:
                    batchfile.write(
                        self._singularity_exec(
                            "blast",
                            f"blastn -db {os.path.dirname(ref)}/{ref_nosuf}  -query {self.finishdir}/assembly/{self.name}_contigs.fasta -out {self.finishdir}/blast_search/{name}/{ref_nosuf}.txt -task megablast -num_threads {self.slurm_header.threads} -outfmt {blast_format}",
                        )
                        + "\n"
                    )
        elif len(file_list) == 1:
            ref_nosuf = re.search(r"(\w+(?:\-\w+)*)\.\w+", os.path.basename(file_list[0])).group(1)
            batchfile.write(
                f"## BLAST {name} search in {self.sample.get('organism').replace('_', ' ').capitalize()}\n"
            )
            batchfile.write(
                self._singularity_exec(
                    "blast",
                    f"blastn -db {os.path.dirname(search_string)}/{ref_nosuf}  -query {self.finishdir}/assembly/{self.name}_contigs.fasta -out {self.finishdir}/blast_search/{name}/{ref_nosuf}.txt -task megablast -num_threads {self.slurm_header.threads} -outfmt {blast_format}",
                )
                + "\n"
            )
        batchfile.write("\n")
        batchfile.close()

    def create_variantsection(self):
        """Creates a job for variant calling based on local alignment"""
        ref = f"{self.folders.genomes}/{self.sample.get('reference')}.fasta"
        localdir = f"{self.finishdir}/alignment"
        outbase = f"{localdir}/{self.name}_{self.sample.get('reference')}"

        # Create run
        batchfile = open(self.batchfile, "a+")
        batchfile.write("# Variant calling based on local alignment\n")
        batchfile.write(f"mkdir {localdir}\n")

        batchfile.write("## Alignment & Deduplication\n")
        batchfile.write(
            self._singularity_exec(
                "bwa",
                f"bwa mem -M -t {self.slurm_header.threads} {ref} {self.concat_files['f']} {self.concat_files['r']} > {outbase}.sam",
            )
            + "\n"
        )
        batchfile.write(
            self._singularity_exec(
                "samtools",
                f"samtools view --threads {self.slurm_header.threads} -b -o {outbase}.bam -T {ref} {outbase}.sam",
            )
            + "\n"
        )
        batchfile.write(
            self._singularity_exec(
                "samtools",
                f"samtools sort --threads {self.slurm_header.threads} -o {outbase}.bam_sort {outbase}.bam",
            )
            + "\n"
        )
        batchfile.write(
            self._singularity_exec(
                "picard",
                f"picard MarkDuplicates I={outbase}.bam_sort O={outbase}.bam_sort_rmdup M={outbase}.stats.dup REMOVE_DUPLICATES=true",
            )
            + "\n"
        )
        batchfile.write(
            self._singularity_exec("samtools", f"samtools index {outbase}.bam_sort_rmdup") + "\n"
        )
        batchfile.write(
            self._singularity_exec("samtools", f"samtools idxstats {outbase}.bam_sort_rmdup")
            + f" &> {outbase}.stats.ref\n"
        )
        # Removal of temp aligment files
        batchfile.write(f"rm {outbase}.bam {outbase}.sam\n")

        batchfile.write("## Primary stats generation\n")
        # Insert stats, dedupped
        batchfile.write(
            self._singularity_exec(
                "picard",
                f"picard CollectInsertSizeMetrics I={outbase}.bam_sort_rmdup O={outbase}.stats.ins H={outbase}.hist.ins",
            )
            + "\n"
        )
        # Coverage
        batchfile.write(
            self._singularity_exec(
                "samtools", f"samtools stats --coverage 1,10000,1 {outbase}.bam_sort_rmdup"
            )
            + f" |grep ^COV | cut -f 2- &> {outbase}.stats.cov\n"
        )
        # Mapped rate, no dedup,dedup in MWGS (trimming has no effect)!
        batchfile.write(
            self._singularity_exec("samtools", f"samtools flagstat {outbase}.bam_sort")
            + f" &> {outbase}.stats.map\n"
        )
        # Total reads, no dedup,dedup in MWGS (trimming has no effect)!
        batchfile.write(
            self._singularity_exec("samtools", f"samtools view -c {outbase}.bam_sort")
            + f" &> {outbase}.stats.raw\n"
        )

        batchfile.write("\n\n")
        batchfile.close()

    def create_preprocsection(self):
        """Concatinates data, possibly trims it, then makes the unstranded reads usable"""
        forward = list()
        reverse = list()

        trimdir = f"{self.finishdir}/trimmed"
        files = self.verify_fastq()
        batchfile = open(self.batchfile, "a+")
        batchfile.write("#Trimmomatic section\n")
        batchfile.write(f"mkdir {trimdir}\n")

        batchfile.write("##Pre-concatination\n")
        for file in files:
            fullfile = f"{self.indir}/{file}"
            # Even indexes = Forward
            if not files.index(file) % 2:
                forward.append(fullfile)
            elif files.index(file) % 2:
                reverse.append(fullfile)
        outfile = files[0].split("_")[0]

        self.concat_files["f"] = f"{self.finishdir}/trimmed/{self.name}_forward_reads.fastq.gz"
        self.concat_files["r"] = f"{self.finishdir}/trimmed/{self.name}_reverse_reads.fastq.gz"
        batchfile.write(f"cat {' '.join(forward)} > {self.concat_files.get('f')}\n")
        batchfile.write(f"cat {' '.join(reverse)} > {self.concat_files.get('r')}\n")

        if self.trimmed:
            fp = f"{trimdir}/{outfile}_trim_front_pair.fastq.gz"
            fu = f"{trimdir}/{outfile}_trim_front_unpair.fastq.gz"
            rp = f"{trimdir}/{outfile}_trim_rev_pair.fastq.gz"
            ru = f"{trimdir}/{outfile}_trim_rev_unpair.fastq.gz"
            batchfile.write("##Trimming section\n")
            batchfile.write(
                self._singularity_exec(
                    "trimmomatic",
                    (
                        f"trimmomatic PE -threads {self.slurm_header.threads}"
                        f" -phred33 {self.concat_files.get('f')} {self.concat_files.get('r')}"
                        f" {fp} {fu} {rp} {ru}"
                        f" ILLUMINACLIP:{self.singularity.trimmomatic_adapters}/NexteraPE-PE.fa:2:30:10"
                        " LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
                    ),
                )
                + "\n"
            )

            batchfile.write("## Interlaced trimmed files\n")
            self.concat_files["f"] = fp
            self.concat_files["r"] = rp
            self.concat_files["i"] = f"{trimdir}/{outfile}_trim_unpair.fastq.gz"

            batchfile.write(f"cat {' '.join([fu, ru])} >> {self.concat_files.get('i')}\n")
        batchfile.write("\n")
        batchfile.close()

    def create_assemblystats_section(self):
        batchfile = open(self.batchfile, "a+")
        batchfile.write("# QUAST QC metrics\n")
        batchfile.write(f"mkdir {self.finishdir}/assembly/quast\n")
        batchfile.write(
            self._singularity_exec(
                "quast",
                f"quast.py {self.finishdir}/assembly/{self.name}_contigs.fasta -o {self.finishdir}/assembly/quast",
            )
            + "\n"
        )
        batchfile.write(
            f"mv {self.finishdir}/assembly/quast/report.tsv {self.finishdir}/assembly/quast/{self.name}_report.tsv\n\n"
        )
        batchfile.close()

    def create_snpsection(self):
        snplist = self.filelist.copy()
        batchfile = open(self.batchfile, "a+")
        name = ""

        # VCFTools filters:
        vcffilter = "--minQ 30 --thin 50 --minDP 3 --min-meanDP 20"
        # BCFTools filters:
        bcffilter = "GL[0]<-500 & GL[1]=0 & QR/RO>30 & QA/AO>30 & QUAL>5000 & ODDS>1100 & GQ>140 & DP>100 & MQM>59 & SAP<15 & PAIRED>0.9 & EPP>3"

        for item in snplist:
            if item.count("/") >= 2:
                name = item.split("/")[-2]
            if "_" in name:
                name = name.split("_")[0]
            batchfile.write(f"# Basecalling for sample {name}\n")
            ref = f"{self.folders.genomes}/{self.sample.get('reference')}.fasta"
            outbase = f"{item}/{name}_{self.sample.get('reference')}"
            batchfile.write(
                f"samtools view -h -q 1 -F 4 -F 256 {outbase}.bam_sort_rmdup | grep -v XA:Z | grep -v SA:Z| samtools view -b - > {self.finishdir}/{name}.unique\n"
            )
            batchfile.write(
                f"freebayes -= --pvar 0.7 -j -J --standard-filters -C 6 --min-coverage 30 --ploidy 1 -f {ref} -b {self.finishdir}/{name}.unique -v {self.finishdir}/{name}.vcf\n"
            )
            batchfile.write(
                f"bcftools view {self.finishdir}/{name}.vcf -o {self.finishdir}/{name}.bcf.gz -O b --exclude-uncalled --types snps\n"
            )
            batchfile.write(f"bcftools index {self.finishdir}/{name}.bcf.gz\n")
            batchfile.write("\n")

            batchfile.write(
                f"vcftools --bcf {self.finishdir}/{name}.bcf.gz {vcffilter} --remove-filtered-all --recode-INFO-all --recode-bcf --out {self.finishdir}/{name}\n"
            )
            batchfile.write(
                f'bcftools view {self.finishdir}/{name}.recode.bcf -i "{bcffilter}" -o {self.finishdir}/{name}.recode.bcf.gz -O b --exclude-uncalled --types snps\n'
            )
            batchfile.write(f"bcftools index {self.finishdir}/{name}.recode.bcf.gz\n\n")

        batchfile.write("# SNP pair-wise distance\n")
        batchfile.write(f"touch {self.finishdir}/stats.out\n")
        while len(snplist) > 1:
            nameOne = ""
            nameTwo = ""
            top = snplist.pop(0)
            if top.count("/") >= 2:
                nameOne = top.split("/")[-2]
            if "_" in nameOne:
                nameOne = nameOne.split("_")[0]
            for entry in snplist:
                if entry.count("/") >= 2:
                    nameTwo = entry.split("/")[-2]
                if "_" in nameTwo:
                    nameTwo = nameTwo.split("_")[0]

                pair = f"{nameOne}_{nameTwo}"
                batchfile.write(
                    f"bcftools isec {self.finishdir}/{nameOne}.recode.bcf.gz {self.finishdir}/{nameTwo}.recode.bcf.gz -n=1 -c all -p {self.finishdir}/tmp -O b\n"
                )
                batchfile.write(
                    f"bcftools merge -O b -o {self.finishdir}/{pair}.bcf.gz --force-samples {self.finishdir}/tmp/0000.bcf {self.finishdir}/tmp/0001.bcf\n"
                )
                batchfile.write(f"bcftools index {self.finishdir}/{pair}.bcf.gz\n")

                batchfile.write(
                    f"echo {pair} $( bcftools stats {self.finishdir}/{pair}.bcf.gz |grep SNPs: | cut -d $'\\t' -f4 ) >> {self.finishdir}/stats.out\n"
                )
                batchfile.write("\n")
        batchfile.close()

    def create_collection(self):
        """Creates collection entry in database"""
        if self.db_pusher.read_exists("Collections", {"ID_collection": self.name}):
            self.db_pusher.delete_records(name=self.name, type="Collections")
            for sample in self.pool:
                self.db_pusher.add_rec(
                    {"ID_collection": self.name, "CG_ID_sample": sample}, "Collections"
                )

        addedprojs = list()
        for sample in self.pool:
            proj = re.search(r"(\w+)A(?:\w+)", sample).group(1)
            if proj not in addedprojs:
                self.create_project(proj)
                addedprojs.append(proj)

    def create_project(self, name):
        """Creates project in database"""
        proj_col = dict()
        proj_col["CG_ID_project"] = name
        proj_col["Customer_ID_project"] = self.sample.get("Customer_ID_project")
        proj_col["Customer_ID"] = self.sample.get("Customer_ID")
        self.db_pusher.add_rec(proj_col, "Projects")
        self.db_pusher.upd_rec({"CG_ID_project": name}, "Projects", proj_col)

    def create_sample(self, name):
        """Creates sample in database"""
        try:
            sample_col = self.db_pusher.read_columns("Samples")
            sample_col["CG_ID_sample"] = self.sample.get("CG_ID_sample")
            sample_col["CG_ID_project"] = self.sample.get("CG_ID_project")
            sample_col["Customer_ID_sample"] = self.sample.get("Customer_ID_sample")
            sample_col["reference_genome"] = self.sample.get("reference")
            sample_col["reference_length"] = self.sample.get("reference_length")
            sample_col["date_analysis"] = self.dt
            sample_col["organism"] = self.sample.get("organism")
            sample_col["application_tag"] = self.sample.get("application_tag")
            sample_col["priority"] = self.sample.get("priority")
            sample_col["date_arrival"] = datetime.strptime(
                self.sample.get("date_arrival"), "%Y-%m-%d %H:%M:%S"
            )
            sample_col["date_sequencing"] = datetime.strptime(
                self.sample.get("date_sequencing"), "%Y-%m-%d %H:%M:%S"
            )
            sample_col["date_libprep"] = datetime.strptime(
                self.sample.get("date_libprep"), "%Y-%m-%d %H:%M:%S"
            )
            sample_col["method_libprep"] = self.sample.get("method_libprep")
            sample_col["method_sequencing"] = self.sample.get("method_sequencing")
            # self.db_pusher.delete_records(sample_col['CG_ID_sample'], 'sample')
            self.db_pusher.add_rec(sample_col, "Samples")
        except Exception:
            self.logger.error(f"Unable to add sample {self.name} to database")

    def project_job(self, single_sample=False):
        if self.dry:
            dry = True
        else:
            dry = False
        jobarray = list()
        if not os.path.exists(self.finishdir):
            os.makedirs(self.finishdir)
        self.create_version_file(finishdir=self.finishdir)
        # Loads project level info.
        try:
            if single_sample:
                self.create_project(self.sample.get("CG_ID_project"))
            elif self.pool:
                self.create_collection()
            else:
                self.create_project(self.name)
        except Exception:
            self.logger.error(f"LIMS interaction failed. Unable to read/write project {self.name}")
        # Writes the job creation sbatch
        if single_sample:
            try:
                self.sample_job()
                headerargs = self.get_headerargs()
                outfile = self.get_sbatch()
                bash_cmd = f"sbatch {headerargs} {outfile}"
                if not dry and outfile != "":
                    samproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                    output, error = samproc.communicate()
                    jobno = re.search(r"(\d+)", str(output)).group(0)
                    jobarray.append(jobno)
                else:
                    self.logger.info(f"Suppressed command: {bash_cmd}")
            except Exception:
                self.logger.error(f"Unable to analyze single sample {self.name}")
        else:
            for ldir in glob.glob(f"{self.indir}/*/"):
                ldir = os.path.basename(os.path.normpath(ldir))
                try:
                    sample_in = f"{self.indir}/{ldir}"
                    sample_out = f"{self.finishdir}/{ldir}"
                    local_sampleinfo = [p for p in self.sampleinfo if p["CG_ID_sample"] == ldir]
                    if local_sampleinfo == []:
                        raise Exception(f"Sample {ldir} has no counterpart in json file")
                    else:
                        local_sampleinfo = local_sampleinfo[0]
                    sample_settings = dict(self.run_settings)
                    sample_settings["input"] = sample_in
                    sample_settings["finishdir"] = sample_out
                    sample_settings["timestamp"] = self.now
                    sample_instance = Job_Creator(
                        log=self.logger,
                        folders=self.folders,
                        slurm_header=self.slurm_header,
                        regex=self.regex,
                        dry=self.dry,
                        config_path=self.config_path,
                        threshold=self.threshold,
                        pubmlst=self.pubmlst,
                        pasteur=self.pasteur,
                        singularity=self.singularity,
                        containers=self.containers,
                        sampleinfo=local_sampleinfo,
                        run_settings=sample_settings,
                    )
                    sample_instance.sample_job()
                    headerargs = sample_instance.get_headerargs()
                    outfile = ""
                    if os.path.isfile(sample_instance.get_sbatch()):
                        outfile = sample_instance.get_sbatch()
                    bash_cmd = f"sbatch {headerargs} {outfile}"
                    if not dry and outfile != "":
                        projproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                        output, error = projproc.communicate()
                        jobno = re.search(r"(\d+)", str(output)).group(0)
                        jobarray.append(jobno)
                    else:
                        self.logger.info(f"Suppressed command: {bash_cmd}")
                except Exception:
                    pass
        if not dry:
            self.finish_job(jobarray, single_sample)

    def _write_mailjob(self, mailfile: str, report: str, custom_conf: str) -> None:
        """Write the mailjob.sh script that runs `microsalt utils finish` after all jobs complete."""
        _ep = next(
            ep for ep in entry_points(group="console_scripts") if ep.value == "microSALT.cli:root"
        )
        microsalt_bin = Path(sys.executable).parent / _ep.name
        with open(mailfile, "w+") as mb:
            mb.write("#!/usr/bin/env bash\n\n")
            mb.write("#Uploading of results to database and production of report\n")
            finish_cmd = (
                f"{microsalt_bin} --config {self.config_path} utils finish {self.finishdir}/sampleinfo.json "
                f"--input {self.finishdir} "
                f"--email {self.regex.mail_recipient} "
                f"--report {report} "
                f"{custom_conf}\n"
            )
            mb.write(finish_cmd)
            mb.write(f"touch {self.finishdir}/run_complete.out\n")

    def finish_job(self, joblist, single_sample=False):
        """Uploads data and sends an email once all analysis jobs are complete."""
        report = "default"
        if self.qc_only:
            report = "qc"
        custom_conf = ""
        if self.config_path:
            custom_conf = f"--config {self.config_path}"

        process = subprocess.Popen("id -un".split(), stdout=subprocess.PIPE)
        user, error = process.communicate()
        user = str(user).replace(".", " ").title()
        # if not os.path.exists(self.finishdir):
        #  os.makedirs(self.finishdir)

        startfile = f"{self.finishdir}/run_started.out"
        configfile = f"{self.finishdir}/config.log"
        mailfile = f"{self.finishdir}/mailjob.sh"
        samplefile = f"{self.finishdir}/sampleinfo.json"
        with open(samplefile, "w+") as outfile:
            json.dump(self.sampleinfo, outfile)
        with open(startfile, "w+") as sb:
            sb.write("#!/usr/bin/env bash\n")
        with open(configfile, "w+") as cb:
            configout = {
                "folders": self.folders.model_dump(),
                "slurm_header": self.slurm_header.model_dump(),
                "regex": self.regex.model_dump(),
                "threshold": self.threshold.model_dump(),
            }
            cb.write(f"ANALYSIS STARTED BY: {user}\n")
            cb.write(json.dumps(configout, indent=2, separators=(",", ":")))

        self._write_mailjob(mailfile, report, custom_conf)

        massagedJobs = list()
        final = ":".join(joblist)
        # Create subtracker if more than 50 samples
        maxlen = 50
        if len(joblist) > maxlen:
            i = 1
            while i <= len(joblist):
                if i + maxlen < len(joblist):
                    massagedJobs.append(":".join(joblist[i - 1 : i + maxlen - 1]))
                else:
                    massagedJobs.append(":".join(joblist[i - 1 : -1]))
                i += maxlen
            for entry in massagedJobs:
                if massagedJobs.index(entry) < len(massagedJobs) - 1:
                    head = f"-A {self.slurm_header.project} -p core -n 1 -t 00:00:10 -J {self.slurm_header.job_prefix}_{self.name}_SUBTRACKER --qos {self.slurm_header.qos} --dependency=afterany:{entry}"
                    bash_cmd = f"sbatch {head} {startfile}"
                    mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                    output, error = mailproc.communicate()
                    jobno = re.search(r"(\d+)", str(output)).group(0)
                    massagedJobs[massagedJobs.index(entry) + 1] += f":{jobno}"
                else:
                    final = entry
                    break

        head = (
            f"-A {self.slurm_header.project} -p core -n 1 -t 6:00:00 "
            f"-J {self.slurm_header.job_prefix}_{self.name}_MAILJOB "
            f"--qos {self.slurm_header.qos} --open-mode append "
            f"--dependency=afterany:{final} --output {self.finishdir}/mailjob.out"
        )
        bash_cmd = f"sbatch {head} {mailfile}"
        mailproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
        output, error = mailproc.communicate()

        try:
            jobno = str(re.search(r"(\d+)", str(output)).group(0))
            joblist.append(jobno)
        except Exception:
            self.logger.info(f"Unable to grab SLURMID for {self.name}")

        try:
            # Generates file with all slurm ids
            slurmname = f"{self.name}_slurm_ids.yaml"
            slurmreport_storedir = Path(self.folders.reports, "trailblazer", slurmname)
            slurmreport_workdir = Path(self.finishdir, slurmname)
            with slurmreport_workdir.open("w") as slurmreport_file:
                yaml.safe_dump(
                    data={"jobs": [str(job) for job in joblist]},
                    stream=slurmreport_file,
                )
            self.logger.info(f"Dump to {slurmreport_workdir} successful")

            self.logger.info(f"Copying slurm report file to {slurmreport_storedir}")
            shutil.copyfile(slurmreport_workdir, slurmreport_storedir)
            self.logger.info("Copy successful.")
        except Exception as e:
            self.logger.info("Unable to generate Trailblazer slurm report file")
            self.logger.error(f"Error while generating Trailblazer slurm report file:\n {e}")

    def sample_job(self):
        """Writes necessary sbatch job for each individual sample"""
        try:
            if not os.path.exists(self.finishdir):
                os.makedirs(self.finishdir)
            try:
                # This is one job
                self.batchfile = f"{self.finishdir}/runfile.sbatch"
                batchfile = open(self.batchfile, "w+")
                batchfile.write("#!/bin/sh\n\n")
                batchfile.write(f"mkdir -p {self.finishdir}\n")
                batchfile.close()
                self.create_preprocsection()
                self.create_variantsection()
                if not self.qc_only:
                    self.create_assemblysection()
                    self.create_assemblystats_section()
                    self.create_blast_search()
                batchfile = open(self.batchfile, "a+")
                batchfile.close()

                self.logger.info(
                    f"Created runfile for sample {self.name} in folder {self.finishdir}"
                )
            except Exception:
                raise
            try:
                self.create_sample(self.name)
            except Exception:
                self.logger.error(f"Unable to access LIMS info for sample {self.name}")
        except Exception as e:
            self.logger.error(f"Unable to create job for sample {self.name}\nSource: {e!s}")
            shutil.rmtree(self.finishdir, ignore_errors=True)
            raise

    def create_blast_search(self):
        reforganism = self.ref_resolver.organism2reference(self.sample.get("organism"))
        self.batchfile = f"{self.finishdir}/runfile.sbatch"
        batchfile = open(self.batchfile, "a+")
        batchfile.write(f"mkdir -p {self.finishdir}/blast_search\n")
        batchfile.close()
        self.blast_subset(
            "mlst",
            f"{self.folders.references}/{reforganism}/*.tfa",
        )
        self.blast_subset("resistance", f"{self.folders.resistances}/*.fsa")
        if reforganism == "escherichia_coli":
            ss = f"{os.path.dirname(self.folders.expec)}/*{os.path.splitext(self.folders.expec)[1]}"
            self.blast_subset("expec", ss)

    def snp_job(self):
        """Writes a SNP calling job for a set of samples"""
        if not os.path.exists(self.finishdir):
            os.makedirs(self.finishdir)

        self.batchfile = f"{self.finishdir}/runfile.sbatch"
        batchfile = open(self.batchfile, "w+")
        batchfile.write("#!/usr/bin/env bash\n")
        batchfile.write(f"mkdir -p {self.finishdir}\n\n")
        batchfile.close()

        self.create_snpsection()
        batchfile = open(self.batchfile, "a+")
        batchfile.close()

        headerline = f"-A {self.slurm_header.project} -p {self.slurm_header.type} -n 1 -t 24:00:00 -J {self.slurm_header.job_prefix}_{self.name} --qos {self.slurm_header.qos} --output {self.finishdir}/slurm_{self.name}.log"
        outfile = self.get_sbatch()
        bash_cmd = f"sbatch {headerline} {outfile}"
        samproc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
        output, error = samproc.communicate()
