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
            self.logger.error(
                "Unable to open profile folder {}".format(self.config["folders"]["profiles"])
            )

    def add_table(self, file):
        try:
            with open("{}/{}".format(self.config["folders"]["profiles"], file), "r") as fh:
                # Sets profile_* headers
                head = fh.readline()
                head = head.rstrip().split("\t")[:8]  # Only consider the first 8 elements
                index = 0

                header = "Table('profile_{}'.format(file), self.metadata,".format(file)
                while index < len(head):
                    # Set ST as PK
                    if head[index] == "ST":
                        header += "Column(head[{}], SmallInteger, primary_key=True),".format(index)
                    header += "Column(head[{}], SmallInteger),".format(index)
                    index = index + 1
                header += ")"
                p = eval(header)
                self.tables[file] = p
        except Exception as e:
            self.logger.error("Unable to open profile file {}".format(file))


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
            self.logger.error(
                "Unable to open profile folder {}".format(self.config["folders"]["profiles"])
            )

    def add_table(self, file):
        try:
            with open("{}/{}".format(self.config["folders"]["profiles"], file), "r") as fh:
                # Sets profile_* headers
                head = fh.readline()
                head = head.rstrip().split("\t")[:8]  # Only consider the first 8 elements
                index = 0

                header = "Table('novel_{}'.format(file), self.metadata,".format(file)
                while index < len(head):
                    # Set ST as PK
                    if head[index] == "ST":
                        header += "Column(head[{}], SmallInteger, primary_key=True),".format(index)
                    # Set Clonal complex as string
                    header += "Column(head[{}], SmallInteger),".format(index)
                    index = index + 1
                header += ")"
                p = eval(header)
                self.tables[file] = p
        except Exception as e:
            self.logger.error("Unable to open profile file {}".format(file))
