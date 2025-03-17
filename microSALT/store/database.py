from sqlalchemy import SingletonThreadPool, create_engine, MetaData
from sqlalchemy.engine.base import Engine
from sqlalchemy.orm import Session, scoped_session, sessionmaker

ENGINE: Engine = None
SESSION: scoped_session = None


def initialize_database(config: dict):
    global ENGINE, SESSION
    ENGINE = create_engine(
        config["database"]["SQLALCHEMY_DATABASE_URI"], poolclass=SingletonThreadPool, connect_args={"check_same_thread": False, "timeout": 15}
    )
    session_factory = sessionmaker(bind=ENGINE)
    SESSION = scoped_session(session_factory)


def get_engine() -> Engine:
    if not ENGINE:
        raise ValueError("Database engine not initialized")
    return ENGINE


def get_session() -> Session:
    if not SESSION:
        raise ValueError("Database session not initialized")
    return SESSION


def drop_all_tables() -> None:
    """Create all tables in status db."""
    engine = get_engine()
    metadata = MetaData()
    metadata.reflect(engine)
    metadata.drop_all(engine)
    print("All tables dropped")

