import logging
import os

__version__ = "4.3.0"

logger = None

logging_levels = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}


def setup_logger(logging_level: str, log_file: str) -> None:
    global logger
    if logging_level not in logging_levels:
        raise ValueError(
            f"Invalid logging level: {logging_level}. Choose from {list(logging_levels.keys())}."
        )
    logger = logging.getLogger("main_logger")
    logger.setLevel(logging_levels[logging_level])
    ch = logging.StreamHandler()
    ch.setLevel(logging_levels[logging_level])
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S"
    )
    ch.setFormatter(formatter)
    fh = logging.FileHandler(os.path.expanduser(log_file))
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
