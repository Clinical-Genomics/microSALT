from sqlalchemy import Engine, create_engine
from sqlalchemy.orm import Session, scoped_session, sessionmaker

from microSALT.exc.exceptions import MicroSALTError

session: scoped_session | None = None
engine: Engine | None = None


def initialize_database(db_uri: str) -> None:
    """Initialize the SQLAlchemy engine and session for status db."""
    global ENGINE, SESSION

    engine = create_engine(db_uri, pool_pre_ping=True)
    session_factory = sessionmaker(engine)
    session = scoped_session(session_factory)


def get_session() -> Session:
    """Get a SQLAlchemy session with a connection to status db."""
    if not session:
        raise MicroSALTError("Database not initialised")
    return session


def get_scoped_session_registry() -> scoped_session | None:
    """Get the scoped session registry for status db."""
    return session


def get_engine() -> Engine:
    """Get the SQLAlchemy engine with a connection to status db."""
    if not engine:
        raise MicroSALTError("Database not initialised")
    return engine
