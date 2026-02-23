#!/usr/bin/env python
# data_files are not supported in pyproject.toml with setuptools; kept here.
from setuptools import setup

setup(
    data_files=[
        ("expec", ["unique_references/ExPEC.fsa"]),
        ("logos", ["artwork/microsalt.jpg", "artwork/swedac.jpg"]),
        (
            "testdata",
            [
                "tests/testdata/sampleinfo_samples.json",
                "tests/testdata/sampleinfo_mlst.json",
                "tests/testdata/sampleinfo_projects.json",
                "tests/testdata/sampleinfo_reports.json",
                "tests/testdata/sampleinfo_expec.json",
                "tests/testdata/sampleinfo_resistance.json",
                "tests/testdata/quast_results.tsv",
                "tests/testdata/blast_single_resistance.txt",
                "tests/testdata/blast_single_loci.txt",
                "tests/testdata/alignment.stats.ref",
                "tests/testdata/alignment.stats.raw",
                "tests/testdata/alignment.stats.map",
                "tests/testdata/alignment.stats.ins",
                "tests/testdata/alignment.stats.dup",
                "tests/testdata/alignment.stats.cov",
                "configExample.json",
            ],
        ),
        (
            "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A1/",
            [
                "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A1/dummysequences_1.fastq.gz",
                "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A1/dummysequences_2.fastq.gz",
            ],
        ),
        (
            "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A2/",
            [
                "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A2/dummysequences_1.fastq.gz",
                "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A2/dummysequences_2.fastq.gz",
            ],
        ),
        (
            "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A3/",
            [
                "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A3/dummysequences_1.fastq.gz",
                "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A3/dummysequences_2.fastq.gz",
            ],
        ),
        (
            "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A4/",
            [
                "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A4/dummysequences_1.fastq.gz",
                "tests/testdata/AAA1234_2000.1.2_3.4.5/AAA1234A4/dummysequences_2.fastq.gz",
            ],
        ),
    ],
)
