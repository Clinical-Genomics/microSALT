"""Table definitions for profile databases.

Profile tables (profile_* and novel_*) cannot use the declarative ORM because
their names and column sets are determined at runtime by scanning the profiles
folder on disk.  SQLAlchemy Core Table objects are used instead.

By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import os
from sqlalchemy import Column, SmallInteger, Table


class ProfileTable:
    """Builds a dict of SQLAlchemy Core Table objects from the profiles folder.

    Each file in the folder produces one table whose name is
    ``{prefix}{filename}``.  The first eight tab-separated fields of the
    file header become the column names; the column named ``ST`` is used as
    the primary key.

    Args:
        prefix: Table name prefix, e.g. ``"profile_"`` or ``"novel_"``.
        metadata: The shared SQLAlchemy MetaData instance.
        config: Application config dict (must contain ``folders.profiles``).
        log: Logger instance.
    """

    def __init__(self, prefix: str, metadata, config, log):
        self.tables = dict()
        self.prefix = prefix
        self.metadata = metadata
        self.config = config
        self.logger = log
        try:
            for filename in os.listdir(self.config["folders"]["profiles"]):
                self._add_table(filename)
        except Exception:
            self.logger.error(
                f"Unable to open profile folder {self.config['folders']['profiles']}"
            )

    def _add_table(self, filename: str) -> None:
        try:
            with open(f"{self.config['folders']['profiles']}/{filename}", "r") as fh:
                head = fh.readline().rstrip().split("\t")[:8]
            columns = [
                Column(col, SmallInteger, primary_key=(col == "ST"))
                for col in head
            ]
            table = Table(f"{self.prefix}{filename}", self.metadata, *columns)
            self.tables[filename] = table
        except Exception:
            self.logger.error(f"Unable to open profile file {filename}")
