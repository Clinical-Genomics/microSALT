#!/usr/bin/env python
"""Migrate ORM table data from an old SQLite database to a new MySQL database.

Tables migrated (insertion order respects FK constraints):
  projects → samples → seq_types → resistances → expacs
  → reports → collections → versions → system_locks

Profile/Novel tables (managed outside of the ORM declarative base) are
intentionally excluded.

Usage:
    python scripts/migrate_sqlite_to_mysql.py \\
        --sqlite  sqlite:////path/to/old/microsalt.db \\
        --mysql   mysql+pymysql://user:pass@host/dbname

The script is idempotent: rows already present in the target (matched by
primary key) are skipped rather than duplicated or overwritten.
"""

import argparse
import sys

from sqlalchemy import String, create_engine, inspect, text
from sqlalchemy.orm import Session

from microSALT.store.orm_models import (
    Base,
    Collections,
    Expacs,
    Projects,
    Reports,
    Resistances,
    Samples,
    Seq_types,
    SystemLock,
    Versions,
)

# Tables in the order they must be inserted to satisfy FK constraints.
TABLES: list[type] = [
    Projects,
    Samples,
    Seq_types,
    Resistances,
    Expacs,
    Reports,
    Collections,
    Versions,
    SystemLock,
]


def _widen_varchar_columns(dst_engine) -> None:
    """ALTER any MySQL VARCHAR column that is narrower than the ORM definition.

    Called after create_all so that tables are guaranteed to exist.
    Handles the case where the schema was already created with an older,
    narrower column definition.
    """
    dst_inspector = inspect(dst_engine)
    existing_tables = dst_inspector.get_table_names()

    with dst_engine.connect() as conn:
        for model in TABLES:
            table_name = model.__tablename__
            if table_name not in existing_tables:
                continue

            actual_cols = {c["name"]: c for c in dst_inspector.get_columns(table_name)}

            for col in inspect(model).mapper.columns:
                if not isinstance(col.type, String):
                    continue
                orm_len = col.type.length
                if orm_len is None:
                    continue
                actual = actual_cols.get(col.name)
                if actual is None:
                    continue
                actual_len = getattr(actual["type"], "length", None)
                if actual_len is not None and actual_len < orm_len:
                    print(
                        f"  Widening {table_name}.{col.name}: "
                        f"VARCHAR({actual_len}) → VARCHAR({orm_len})"
                    )
                    conn.execute(
                        text(
                            f"ALTER TABLE `{table_name}` MODIFY COLUMN"
                            f" `{col.name}` VARCHAR({orm_len})"
                        )
                    )
        conn.commit()


def _columns(model: type) -> list[str]:
    """Return the list of column attribute names for an ORM model."""
    return [c.key for c in inspect(model).mapper.column_attrs]


def _pk_columns(model: type) -> list[str]:
    """Return the primary key column names for an ORM model."""
    return [col.name for col in inspect(model).mapper.primary_key]


def migrate_table(src: Session, dst: Session, model: type) -> tuple[int, int]:
    """Copy rows from src to dst for the given model.

    Returns (inserted, skipped) counts.
    """
    table_name = model.__tablename__
    cols = _columns(model)
    pk_cols = _pk_columns(model)

    inserted = 0
    skipped = 0

    rows = src.query(model).all()
    for row in rows:
        # Check if the row already exists in the destination by PK lookup
        pk_filter = {col: getattr(row, col) for col in pk_cols}
        exists = dst.get(
            model,
            tuple(pk_filter[c] for c in pk_cols)
            if len(pk_cols) > 1
            else next(iter(pk_filter.values())),
        )
        if exists is not None:
            skipped += 1
            continue

        # Build a fresh detached copy so we don't accidentally modify the
        # source session's identity map.
        new_obj = model(**{col: getattr(row, col) for col in cols})
        dst.add(new_obj)
        inserted += 1

    dst.flush()
    return inserted, skipped


def main() -> int:
    parser = argparse.ArgumentParser(description="Migrate microSALT ORM data from SQLite to MySQL.")
    parser.add_argument(
        "--sqlite",
        required=True,
        metavar="PATH_OR_URI",
        help="Path to the source SQLite file, or a full SQLAlchemy URI  (e.g. sqlite:////path/to/microsalt.db)",
    )
    parser.add_argument(
        "--mysql",
        required=True,
        metavar="URI",
        help="SQLAlchemy URI for the target MySQL DB  (e.g. mysql+pymysql://user:pass@host/db)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Read source data and report counts without writing to MySQL.",
    )
    args = parser.parse_args()

    # Accept a bare file path as well as a full SQLAlchemy URI
    sqlite_uri = args.sqlite
    if not sqlite_uri.startswith("sqlite:"):
        import os

        sqlite_uri = "sqlite:///" + os.path.abspath(sqlite_uri)

    print(f"Source : {sqlite_uri}")
    print(f"Target : {args.mysql}")
    if args.dry_run:
        print("DRY RUN — no data will be written.\n")

    src_engine = create_engine(sqlite_uri, pool_pre_ping=True)
    dst_engine = create_engine(args.mysql, pool_pre_ping=True)

    # Ensure all ORM tables exist in the destination, then widen any columns
    # that are narrower in MySQL than the current ORM definition.
    if not args.dry_run:
        Base.metadata.create_all(dst_engine)
        _widen_varchar_columns(dst_engine)

    total_inserted = 0
    total_skipped = 0

    with Session(src_engine) as src_session, Session(dst_engine) as dst_session:
        src_inspector = inspect(src_engine)
        for model in TABLES:
            table_name = model.__tablename__

            # Check whether the table exists in the source at all
            if table_name not in src_inspector.get_table_names():
                print(f"  {table_name:<20} — not present in source, skipping")
                continue

            row_count = src_session.query(model).count()
            if row_count == 0:
                print(f"  {table_name:<20} — 0 rows in source, skipping")
                continue

            if args.dry_run:
                print(f"  {table_name:<20} — {row_count} rows would be processed")
                total_inserted += row_count
                continue

            try:
                inserted, skipped = migrate_table(src_session, dst_session, model)
                print(
                    f"  {table_name:<20} — inserted {inserted}, skipped {skipped} (already present)"
                )
                total_inserted += inserted
                total_skipped += skipped
            except Exception as exc:
                dst_session.rollback()
                print(f"  {table_name:<20} — ERROR: {exc}", file=sys.stderr)
                print("Rolling back entire transaction.", file=sys.stderr)
                return 1

        if not args.dry_run:
            dst_session.commit()

    print(f"\nDone. Inserted {total_inserted} rows, skipped {total_skipped} duplicates.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
