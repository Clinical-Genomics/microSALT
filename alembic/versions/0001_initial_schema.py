"""initial schema

Revision ID: 20260424_0001
Revises:
Create Date: 2026-04-24

"""

from typing import Sequence, Union

from alembic import op

# revision identifiers, used by Alembic.
revision: str = "20260424_0001"
down_revision: Union[str, Sequence[str], None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    from microSALT.store.orm_models import Base

    Base.metadata.create_all(bind=op.get_bind())


def downgrade() -> None:
    from microSALT.store.orm_models import Base

    Base.metadata.drop_all(bind=op.get_bind())
