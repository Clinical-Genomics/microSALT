#!/usr/bin/env python
from microSALT import __version__
from setuptools import setup, find_packages

version = __version__

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []

setup(
    name="microSALT",
    version=version,
    long_description=__doc__,
    url="https://github.com/Clinical-Genomics/microSALT",
    author="Isak Sylvin",
    author_email='isak.sylvin@scilifelab.se',
    install_requires=install_requires,
    packages=find_packages(),
    include_package_data=True,
    zip_safe=False,
    data_files=[('expec', ['unique_references/ExPEC.fsa']),
                ('logos', ['artwork/microsalt.jpg', 'artwork/swedac.jpg']),
                ('testdata', ['tests/testdata/sampleinfo_samples.json','tests/testdata/sampleinfo_mlst.json','tests/testdata/quast_results.tsv', 'tests/testdata/blast_single_resistance.txt','tests/testdata/blast_single_loci.txt','tests/testdata/alignment.stats.ref','tests/testdata/alignment.stats.raw','tests/testdata/alignment.stats.map','tests/testdata/alignment.stats.ins','tests/testdata/alignment.stats.dup','tests/testdata/alignment.stats.cov','configExample.json']) 
],
    entry_points={
        'console_scripts': ['microSALT=microSALT.cli:root'],
    },
)

