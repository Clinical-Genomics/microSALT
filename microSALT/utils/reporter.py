"""Generates various reports by tapping into Flask and mySQL
By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import json
import os
import smtplib
import socket
import sys
import time
from datetime import datetime
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from shutil import copyfile

import yaml

from microSALT.server.views import (
    alignment_page,
    gen_collectiondata,
    gen_reportdata,
    STtracker_page,
    typing_page,
)
from microSALT.store.db_manipulator import DB_Manipulator


class Reporter:
    def __init__(self, config, log, sampleinfo={}, name="", output="", collection=False):
        self.db_pusher = DB_Manipulator(config, log)
        self.name = name
        self.collection = collection
        if output == "":
            self.output = os.getcwd()
        else:
            self.output = output + "/"
        self.config = config
        self.logger = log
        self.attachments = list()
        self.filedict = dict()
        self.error = False
        self.dt = datetime.now()
        self.now = time.strftime(
            f"{self.dt.year}.{self.dt.month}.{self.dt.day}_{self.dt.hour}.{self.dt.minute}.{self.dt.second}"
        )

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
            self.name = self.sampleinfo.get("CG_ID_project")
            self.sample = self.sampleinfo

    def create_subfolders(self):
        os.makedirs(f"{self.config['folders']['reports']}/deliverables", exist_ok=True)
        os.makedirs(f"{self.config['folders']['reports']}/json", exist_ok=True)
        os.makedirs(f"{self.config['folders']['reports']}/analysis", exist_ok=True)

    def report(self, type="default", customer="all"):
        self.create_subfolders()
        if type in ["default", "typing", "qc"]:
            # Only typing and qc reports are version controlled
            self.gen_version(self.name)
        if type in ["default", "typing", "qc", "st_update"]:
            if type == "default":
                self.gen_typing()
                self.gen_qc()
                self.gen_json(silent=True)
                self.gen_delivery()
            elif type == "typing":
                self.gen_typing()
            elif type == "qc":
                self.gen_qc()
            elif type == "st_update":
                self.gen_STtracker(customer)
        elif type in ["json_dump", "motif_overview"]:
            if type == "json_dump":
                self.gen_json()
                self.gen_delivery()
            elif type == "motif_overview":
                self.gen_motif(motif="resistance")
                self.gen_motif(motif="expec")
        else:
            raise Exception("Report function recieved invalid format")
        self.mail()
        # If no output dir is specified; Don't store report locally. Rely on e-mail
        if not self.output == "" or self.output == os.getcwd():
            for k, v in self.filedict.items():
                if v == "":
                    os.remove(k)
                else:
                    copyfile(k, v)

    def gen_version(self, name):
        self.db_pusher.get_report(name)
        self.db_pusher.set_report(name)

    def gen_STtracker(self, customer="all", silent=False):
        self.name = "Sequence Type Update"
        try:
            content = STtracker_page(customer)
            outname = f"{self.output}/ST_updates_{self.now}.html"
            outfile = open(outname, "wb")
            outfile.write(content.encode("utf8"))
            outfile.close()
            self.filedict[outname] = ""
            if not silent:
                self.attachments.append(outname)
        except Exception:
            self.logger.error("Failed to generate ST tracker report")
            self.error = True

    def gen_qc(self, silent=False):
        try:
            last_version = self.db_pusher.get_report(self.name).version
        except Exception:
            self.logger.error(f"Project {self.name} does not exist")
            sys.exit(-1)
        try:
            content = alignment_page(self.name)
            outfile = f"{self.sample.get('Customer_ID_project')}_QC_{last_version}.html"
            local = f"{self.output}/{outfile}"
            output = f"{self.config['folders']['reports']}/analysis/{outfile}"

            with open(output, "wb") as f:
                f.write(content.encode("utf8"))

            if os.path.isfile(output):
                self.filedict[output] = local
                if not silent:
                    self.attachments.append(output)
        except Exception:
            self.logger.error("Failed to generate QC report")
            self.error = True

    def gen_typing(self, silent=False):
        try:
            last_version = self.db_pusher.get_report(self.name).version
        except Exception:
            self.logger.error(f"Project {self.name} does not exist")
            sys.exit(-1)
        try:
            content = typing_page(self.name, "all")
            outfile = f"{self.sample.get('Customer_ID_project')}_Typing_{last_version}.html"
            local = f"{self.output}/{outfile}"
            output = f"{self.config['folders']['reports']}/analysis/{outfile}"

            with open(output, "wb") as f:
                f.write(content.encode("utf8"))

            if os.path.isfile(output):
                self.filedict[output] = local
                if not silent:
                    self.attachments.append(output)
        except Exception:
            self.logger.error("Failed to generate typing report")
            self.error = True

    def gen_motif(self, motif="resistance", silent=False):
        if motif not in ["resistance", "expec"]:
            self.logger.error("Invalid motif type specified for gen_motif function")
        if self.collection:
            sample_info = gen_collectiondata(self.name)
        else:
            sample_info = gen_reportdata(self.name)
        output = f"{self.output}/{self.name}_{motif}_{self.now}.csv"

        # Load motif & gene names into dict
        motifdict = dict()
        for s in sample_info["samples"]:
            if motif == "resistance":
                for r in s.resistances:
                    if r.resistance not in motifdict.keys() and r.threshold == "Passed":
                        if r.resistance is None:
                            r.resistance = "None"
                        motifdict[r.resistance] = list()
                    if r.threshold == "Passed" and r.gene not in motifdict[r.resistance]:
                        motifdict[r.resistance].append(r.gene)
            elif motif == "expec":
                for e in s.expacs:
                    if e.virulence not in motifdict.keys() and e.threshold == "Passed":
                        if e.virulence is None:
                            e.virulence = "None"
                        motifdict[e.virulence] = list()
                    if e.threshold == "Passed" and e.gene not in motifdict[e.virulence]:
                        motifdict[e.virulence].append(e.gene)
        for k, v in motifdict.items():
            motifdict[k] = sorted(v)

        # Top 2 Header
        sepfix = "sep=,"
        topline = f"Identity {self.config['threshold']['motif_id']}% & Span {self.config['threshold']['motif_span']}%,,,"
        botline = "CG Sample ID,Sample ID,Organism,Sequence Type,Thresholds"
        for k in sorted(motifdict.keys()):
            genes = [""] * len(motifdict[k])
            active_gene = k.replace(",", " &")
            if active_gene == "":
                active_gene = "Uncategorized hits"
            geneholder = ",".join(genes)
            topline += f",,{active_gene}{geneholder}"
            resnames = ",".join(sorted(motifdict[k]))
            botline += f",,{resnames}"

        try:
            excel = open(output, "w+")
            excel.write(f"{sepfix}\n")
            excel.write(f"{topline}\n")
            excel.write(f"{botline}\n")

            # Create each individual row past the 2nd, per iteration
            for s in sample_info["samples"]:
                rowdict = dict()
                pref = f"{s.CG_ID_sample},{s.Customer_ID_sample},{s.organism},{s.ST_status.replace(',', ';')},{s.threshold}"
                # Load single sample
                if motif == "resistance":
                    for r in s.resistances:
                        if r.resistance not in rowdict.keys() and r.threshold == "Passed":
                            rowdict[r.resistance] = dict()
                        if r.threshold == "Passed" and r.gene not in rowdict[r.resistance]:
                            rowdict[r.resistance][r.gene] = r.identity
                elif motif == "expec":
                    for e in s.expacs:
                        if e.virulence not in rowdict.keys() and e.threshold == "Passed":
                            rowdict[e.virulence] = dict()
                        if e.threshold == "Passed" and e.gene not in rowdict[e.virulence]:
                            rowdict[e.virulence][e.gene] = e.identity
                # Compare single sample to all
                hits = ""
                for res in sorted(motifdict.keys()):
                    if res in rowdict.keys():
                        hits += ",1"
                        for gen in sorted(motifdict[res]):
                            hits += ","
                            if gen in rowdict[res].keys():
                                # UPD: Change this to identity of hit
                                hits += f"{rowdict[res][gen]}"
                            else:
                                hits += "0"
                    else:
                        # Commas eq to res + gen length
                        hits += ",0,0"
                        pad = ["0"] * len(motifdict[res])
                        hits += ",".join(pad)

                excel.write(f"{pref}{hits}\n")

            excel.close()
            if os.path.isfile(output):
                self.filedict[output] = ""
                if not silent:
                    self.attachments.append(output)
        except FileNotFoundError as e:
            self.logger.error(
                f"Gen_motif unable to produce excel file. Path {os.path.basename(output)} does not exist"
            )

    def gen_delivery(self):
        deliv = dict()
        deliv["files"] = list()
        last_version = self.db_pusher.get_report(self.name).version
        output = f"{self.config['folders']['reports']}/deliverables/{self.sample.get('Customer_ID_project')}_deliverables.yaml"
        local = f"{self.output}/{self.sample.get('Customer_ID_project')}_deliverables.yaml"

        # Project-wide
        # Sampleinfo
        deliv["files"].append(
            {
                "format": "json",
                "id": str(self.sample.get("Customer_ID_project")),
                "path": f"{self.output}/sampleinfo.json",
                "path_index": "~",
                "step": "analysis",
                "tag": "sampleinfo",
            }
        )
        # QC report
        deliv["files"].append(
            {
                "format": "html",
                "id": str(self.sample.get("Customer_ID_project")),
                "path": f"{self.output}/{self.sample.get('Customer_ID_project')}_QC_{last_version}.html",
                "path_index": "~",
                "step": "result_aggregation",
                "tag": "microsalt-qc",
            }
        )
        # Typing report
        deliv["files"].append(
            {
                "format": "html",
                "id": str(self.sample.get("Customer_ID_project")),
                "path": f"{self.output}/{self.sample.get('Customer_ID_project')}_Typing_{last_version}.html",
                "path_index": "~",
                "step": "result_aggregation",
                "tag": "microsalt-type",
            }
        )
        # Json (vogue) report
        deliv["files"].append(
            {
                "format": "json",
                "id": str(self.sample.get("Customer_ID_project")),
                "path": f"{self.output}/{self.sample.get('CG_ID_project')}.json",
                "path_index": "~",
                "step": "result_aggregation",
                "tag": "microsalt-json",
            }
        )
        # Settings dump
        deliv["files"].append(
            {
                "format": "txt",
                "id": str(self.sample.get("Customer_ID_project")),
                "path": f"{self.output}/config.log",
                "path_index": "~",
                "step": "analysis",
                "tag": "runtime-settings",
            }
        )

        # Version file
        deliv["files"].append(
            {
                "format": "txt",
                "id": str(self.sample.get("Customer_ID_project")),
                "path": f"{self.output}/version.txt",
                "path_index": "~",
                "step": "result_aggregation",
                "tag": "microsalt-version",
            }
        )

        # Sample-wide
        # Single sample
        if self.sampleinfo == self.sample:
            hklist = list()
            hklist.append(self.sampleinfo)
            resultsdir = self.output
        # Project
        else:
            hklist = self.sampleinfo

        for s in hklist:
            if len(hklist) > 1:
                resultsdir = os.path.join(self.output, s["CG_ID_sample"])
            # Contig/Assembly file
            deliv["files"].append(
                {
                    "format": "fasta",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/assembly/{s['CG_ID_sample']}_trimmed_contigs.fasta",
                    "path_index": "~",
                    "step": "assembly",
                    "tag": "assembly",
                }
            )
            # Concat trimmed reads forwards
            deliv["files"].append(
                {
                    "format": "fastq",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/trimmed/{s['CG_ID_sample']}_trim_front_pair.fastq.gz",
                    "path_index": "~",
                    "step": "concatination",
                    "tag": "trimmed-forward-reads",
                }
            )
            # Concat trimmed reads reverse
            deliv["files"].append(
                {
                    "format": "fastq",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/trimmed/{s['CG_ID_sample']}_trim_rev_pair.fastq.gz",
                    "path_index": "~",
                    "step": "concatination",
                    "tag": "trimmed-reverse-reads",
                }
            )
            # Concat trimmed reads unpaired
            deliv["files"].append(
                {
                    "format": "fastq",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/trimmed/{s['CG_ID_sample']}_trim_unpair.fastq.gz",
                    "path_index": "~",
                    "step": "concatination",
                    "tag": "trimmed-unpaired-reads",
                }
            )
            # Slurm dump
            deliv["files"].append(
                {
                    "format": "txt",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/slurm_{s['CG_ID_sample']}.log",
                    "path_index": "~",
                    "step": "analysis",
                    "tag": "logfile",
                }
            )
            # Quast (assembly) qc report
            deliv["files"].append(
                {
                    "format": "tsv",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/assembly/quast/{s['CG_ID_sample']}_report.tsv",
                    "path_index": "~",
                    "step": "assembly",
                    "tag": "quast-results",
                }
            )
            # Alignment (bam, sorted)
            deliv["files"].append(
                {
                    "format": "bam",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/alignment/{s['CG_ID_sample']}_{s['reference']}.bam_sort",
                    "path_index": "~",
                    "step": "alignment",
                    "tag": "reference-alignment-sorted",
                }
            )
            # Alignment (bam, sorted, deduplicated)
            deliv["files"].append(
                {
                    "format": "bam",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/alignment/{s['CG_ID_sample']}_{s['reference']}.bam_sort_rmdup",
                    "path_index": "~",
                    "step": "alignment",
                    "tag": "reference-alignment-deduplicated",
                }
            )
            # Picard insert size stats
            deliv["files"].append(
                {
                    "format": "meta",
                    "id": s["CG_ID_sample"],
                    "path": f"{resultsdir}/alignment/{s['CG_ID_sample']}_{s['reference']}.stats.ins",
                    "path_index": "~",
                    "step": "insertsize_calc",
                    "tag": "picard-insertsize",
                }
            )

        with open(output, "w") as delivfile:
            documents = yaml.dump(deliv, delivfile)

        with open(output, "r") as delivfile:
            postfix = delivfile.read()
        postfix = postfix.replace("'~'", "~")

        with open(output, "w") as delivfile:
            delivfile.write(postfix)

        if os.path.isfile(output):
            self.filedict[output] = local

    def gen_json(self, silent=False):
        report = dict()
        local = f"{self.output}/{self.name}.json"
        output = f"{self.config['folders']['reports']}/json/{self.name}.json"

        sample_info = gen_reportdata(self.name)
        analyses = [
            "blast_pubmlst",
            "quast_assembly",
            "blast_resfinder_resistence",
            "picard_markduplicate",
            "microsalt_samtools_stats",
        ]
        for s in sample_info["samples"]:
            t = dict()

            # Since some apps are too basic to filter irrelevant non-standard values..
            t["ST_status"] = "" if s.ST_status is None or s.ST_status != str(s.ST) else s.ST_status
            t["threshold"] = (
                ""
                if s.threshold is None or s.threshold not in ["Passed", "Failed"]
                else s.threshold
            )
            t["genome_length"] = (
                "" if s.genome_length is None or s.genome_length < 1 else s.genome_length
            )
            t["reference_length"] = (
                "" if s.reference_length is None or s.reference_length < 1 else s.reference_length
            )
            t["gc_percentage"] = (
                "" if s.gc_percentage is None or s.gc_percentage < 0.1 else str(s.gc_percentage)
            )
            t["n50"] = "" if s.n50 is None or s.n50 < 1 else s.n50
            t["contigs"] = "" if s.contigs is None or s.contigs < 1 else s.contigs
            t["insert_size"] = "" if s.insert_size is None or s.insert_size < 1 else s.insert_size
            t["duplication_rate"] = "" if s.duplication_rate is None else s.duplication_rate
            t["total_reads"] = "" if s.total_reads is None or s.total_reads < 1 else s.total_reads
            t["mapped_rate"] = "" if s.mapped_rate is None or s.mapped_rate < 0.1 else s.mapped_rate
            t["average_coverage"] = (
                "" if s.average_coverage is None or s.average_coverage < 0.1 else s.average_coverage
            )
            t["coverage_10x"] = (
                "" if s.coverage_10x is None or s.coverage_10x < 0.1 else s.coverage_10x
            )
            t["coverage_30x"] = (
                "" if s.coverage_30x is None or s.coverage_30x < 0.1 else s.coverage_30x
            )
            t["coverage_50x"] = (
                "" if s.coverage_50x is None or s.coverage_50x < 0.1 else s.coverage_50x
            )
            t["coverage_100x"] = (
                "" if s.coverage_100x is None or s.coverage_100x < 0.1 else s.coverage_100x
            )

            report[s.CG_ID_sample] = dict()
            for a in analyses:
                if a == "blast_resfinder_resistence":
                    report[s.CG_ID_sample][a] = list()
                else:
                    report[s.CG_ID_sample][a] = dict()

            report[s.CG_ID_sample]["blast_pubmlst"] = {
                "sequence_type": t["ST_status"],
                "thresholds": t["threshold"],
            }
            report[s.CG_ID_sample]["quast_assembly"] = {
                "estimated_genome_length": t["genome_length"],
                "gc_percentage": t["gc_percentage"],
                "n50": t["n50"],
                "necessary_contigs": t["contigs"],
            }
            report[s.CG_ID_sample]["picard_markduplicate"] = {
                "insert_size": t["insert_size"],
                "duplication_rate": t["duplication_rate"],
            }
            report[s.CG_ID_sample]["microsalt_samtools_stats"] = {
                "total_reads": t["total_reads"],
                "mapped_rate": t["mapped_rate"],
                "average_coverage": t["average_coverage"],
                "coverage_10x": t["coverage_10x"],
                "coverage_30x": t["coverage_30x"],
                "coverage_50x": t["coverage_50x"],
                "coverage_100x": t["coverage_100x"],
            }

            for r in s.resistances:
                if (
                    r.gene not in report[s.CG_ID_sample]["blast_resfinder_resistence"]
                    and r.threshold == "Passed"
                ):
                    report[s.CG_ID_sample]["blast_resfinder_resistence"].append(r.gene)

        # json.dumps(report) #Dumps the json directly
        try:
            with open(output, "w") as outfile:
                json.dump(report, outfile)

            if os.path.isfile(output):
                self.filedict[output] = local
                if not silent:
                    self.attachments.append(output)
        except FileNotFoundError:
            self.logger.error(
                f"Gen_json unable to produce json file. Path {os.path.basename(output)} does not exist"
            )

    def mail(self):
        msg = MIMEMultipart()
        if not self.error and self.attachments:
            msg["Subject"] = f"{self.name} ({self.attachments[0].split('_')[0]}) Reports"
        else:
            msg["Subject"] = f"{self.name} Failed Generating Report"

        sender = socket.gethostname()
        sender_fixed = f"{os.path.splitext(sender)[0]}.com"
        msg["From"] = sender_fixed

        msg["To"] = self.config["regex"]["mail_recipient"]

        if not self.error:
            for file in self.attachments:
                part = MIMEApplication(open(file).read())
                part.add_header(
                    "Content-Disposition",
                    f'attachment; filename="{os.path.basename(file)}"',
                )
                msg.attach(part)

        s = smtplib.SMTP("localhost")
        s.connect()
        s.sendmail(msg["From"], msg["To"], msg.as_string())
        s.quit()
        self.logger.info(f"Mail containing report sent to {msg['To']} from {msg['From']}")
