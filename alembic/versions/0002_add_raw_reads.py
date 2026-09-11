"""add raw_reads column

Revision ID: 20260903_0002
Revises: 20260424_0001
Create Date: 2026-09-03

"""

from typing import Sequence, Union

import sqlalchemy as sa
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "20260903_0002"
down_revision: Union[str, Sequence[str], None] = "20260424_0001"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    bind = op.get_bind()
    inspector = sa.inspect(bind)
    columns = {col["name"] for col in inspector.get_columns("samples")}
    if "raw_reads" not in columns:
        op.add_column("samples", sa.Column("raw_reads", sa.Integer(), nullable=True))


def downgrade() -> None:
    op.drop_column("samples", "raw_reads")
