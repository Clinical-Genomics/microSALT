import collections
import collections.abc
import json
import logging
import os
import pathlib
import re
import subprocess
import sys
from importlib.resources import files as resource_files

__version__ = "4.3.0"

# Keep track of microSALT installation
wd = os.path.dirname(os.path.realpath(__file__))

# Load configuration
preset_config = ""

logger = None

logging_levels = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}


def setup_logger(logging_level: str, preset_config) -> None:
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
    fh = logging.FileHandler(os.path.expanduser(preset_config["folders"]["log_file"]))
    logger.addHandler(fh)
    logger.addHandler(ch)


default = os.path.join(os.path.dirname(wd), "configExample.json")

if "MICROSALT_CONFIG" in os.environ:
    try:
        envvar = os.environ["MICROSALT_CONFIG"]
        with open(envvar, "r") as conf:
            preset_config = json.load(conf)
    except Exception as e:
        print(f"Config error: {e!s}")
        pass
elif os.path.exists(default):
    try:
        with open(os.path.abspath(default), "r") as conf:
            preset_config = json.load(conf)
    except Exception as e:
        print(f"Config error: {e!s}")
        pass

# Config dependent section:
CONFIG = {}

if preset_config != "":
    try:
        CONFIG = {}

        # Initialize database
        from microSALT.store.database import initialize_database

        initialize_database(preset_config["database"]["SQLALCHEMY_DATABASE_URI"])

        # Add `folders` configuration
        CONFIG["folders"] = preset_config.get("folders", {})

        # Ensure PubMLST configuration is included

        CONFIG["pubmlst"] = preset_config.get("pubmlst", {"client_id": "", "client_secret": ""})

        CONFIG["pasteur"] = preset_config.get("pasteur", {"client_id": "", "client_secret": ""})

        # Add extrapaths to config
        preset_config["folders"]["expec"] = str(
            resource_files("microSALT").joinpath("unique_references", "ExPEC.fsa")
        )
        preset_config["singularity"][
            "trimmomatic_adapters"
        ] = "/opt/conda/share/trimmomatic/adapters/"

        # Initialize logger
        setup_logger(logging_level="INFO", preset_config=preset_config)

        # Create paths mentioned in config
        db_file = re.search(
            "sqlite:///(.+)",
            preset_config["database"]["SQLALCHEMY_DATABASE_URI"],
        ).group(1)
        for entry in preset_config.keys():
            if entry not in ["_comment", "singularity", "genologics"]:
                if isinstance(preset_config[entry], str) and "/" in preset_config[entry]:
                    if not preset_config[entry].startswith("/"):
                        sys.exit(-1)
                    unmade_fldr = os.path.abspath(preset_config[entry])
                    if not pathlib.Path(unmade_fldr).exists():
                        os.makedirs(unmade_fldr)
                        logger.info(f"Created path {unmade_fldr}")

                # level two
                elif isinstance(preset_config[entry], collections.abc.Mapping):
                    for thing in preset_config[entry].keys():
                        if (
                            isinstance(preset_config[entry][thing], str)
                            and "/" in preset_config[entry][thing]
                        ):
                            # Special string, mangling
                            if thing == "log_file":
                                unmade_fldr = os.path.dirname(preset_config[entry][thing])
                                bash_cmd = f"touch {preset_config[entry][thing]}"
                                proc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                                output, error = proc.communicate()
                            elif thing == "SQLALCHEMY_DATABASE_URI":
                                unmade_fldr = os.path.dirname(db_file)
                                bash_cmd = f"touch {db_file}"
                                proc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                                output, error = proc.communicate()
                                if proc.returncode != 0:
                                    logger.error(
                                        "Database writing failed! Invalid user access detected!"
                                    )
                                    sys.exit(-1)
                            else:
                                unmade_fldr = preset_config[entry][thing]
                            if not pathlib.Path(unmade_fldr).exists():
                                os.makedirs(unmade_fldr)
                                logger.info(f"Created path {unmade_fldr}")

    except Exception as e:
        print(f"Config error: {e!s}")
        pass
