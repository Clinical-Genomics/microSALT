"""Table definitions for profiles databases. Bit special since it spawns multiple tables.
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import os
from sqlalchemy import *


class Profiles:
    def __init__(self, metadata, config, log):
        self.tables = dict()
        self.metadata = metadata
        self.config = config
        self.logger = log
        try:
            indata = os.listdir(self.config["folders"]["profiles"])
            for file in indata:
                self.add_table(file)
        except Exception as e:
            self.logger.error(f"Unable to open profile folder {self.config['folders']['profiles']}")

    def add_table(self, file):
        try:
            with open(f"{self.config['folders']['profiles']}/{file}", "r") as fh:
                # Sets profile_* headers
                head = fh.readline()
                head = head.rstrip().split("\t")[:8]  # Only consider the first 8 elements
                columns = []
                for col_name in head:
                    if col_name == "ST":
                        columns.append(Column(col_name, SmallInteger, primary_key=True))
                    else:
                        columns.append(Column(col_name, SmallInteger))
                p = Table(f"profile_{file}", self.metadata, *columns)
                self.tables[file] = p
        except Exception as e:
            self.logger.error(f"Unable to open profile file {file}")


class Novel:
    def __init__(self, metadata, config, log):
        self.tables = dict()
        self.metadata = metadata
        self.config = config
        self.logger = log
        try:
            indata = os.listdir(self.config["folders"]["profiles"])
            for file in indata:
                self.add_table(file)
        except Exception as e:
            self.logger.error(f"Unable to open profile folder {self.config['folders']['profiles']}")

    def add_table(self, file):
        try:
            with open(f"{self.config['folders']['profiles']}/{file}", "r") as fh:
                # Sets profile_* headers
                head = fh.readline()
                head = head.rstrip().split("\t")[:8]  # Only consider the first 8 elements
                columns = []
                for col_name in head:
                    if col_name == "ST":
                        columns.append(Column(col_name, SmallInteger, primary_key=True))
                    # Set Clonal complex as string
                    else:
                        columns.append(Column(col_name, SmallInteger))
                p = Table(f"novel_{file}", self.metadata, *columns)
                self.tables[file] = p
        except Exception as e:
            self.logger.error(f"Unable to open profile file {file}")
