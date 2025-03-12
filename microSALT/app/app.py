import collections
import logging
import json
import os
from pathlib import Path
import re
import subprocess
import sys

from flask import Flask
from distutils.sysconfig import get_python_lib

from microSALT.exc.exc import MissingConfigError


def load_config():
    if "MICROSALT_CONFIG" not in os.environ:
        raise MissingConfigError(
            "No config file found! Please set the environment variable MICROSALT_CONFIG to the path of the config file."
        )
    try:
        with open(os.environ["MICROSALT_CONFIG"], "r") as conf:
            return json.load(conf)
    except Exception as e:
            print(f"Config error: {str(e)}")

    raise MissingConfigError(
        "No config file found! Please set the environment variable MICROSALT_CONFIG to the path of the config file."
    )


def initialize_logger(preset_config):
    logger = logging.getLogger("main_logger")
    logger.setLevel(logging.INFO)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))
    logger.addHandler(ch)

    if "folders" in preset_config and "log_file" in preset_config["folders"]:
        fh = logging.FileHandler(os.path.expanduser(preset_config["folders"]["log_file"]))
        fh.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logger.addHandler(fh)

    return logger


def create_path(path: str, logger: logging.Logger):
    if not Path(path).exists():
        os.makedirs(path)
        logger.info(f"Created path {path}")


def handle_config_entry(entry, value, logger):
    if isinstance(value, str) and "/" in value and entry not in ["genologics"]:
        if not value.startswith("/"):
            sys.exit(-1)
        create_path(os.path.abspath(value), logger)
    elif isinstance(value, collections.Mapping):
        for sub_entry, sub_value in value.items():
            handle_config_entry(sub_entry, sub_value, logger)


def create_paths(preset_config, logger):

    for entry, value in preset_config.items():
        if entry == "_comment":
            continue
        handle_config_entry(entry, value, logger)


def check_database_integrity(db_file, logger):
    cmd = f"sqlite3 {db_file} pragma integrity_check;"
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    if "ok" not in str(output):
        logger.error("Database integrity failed! Lock-state detected!")
        sys.exit(-1)


def _setup_app() -> Flask:
    app = Flask("microSALT", template_folder="server/templates")
    app.config.setdefault("SQLALCHEMY_DATABASE_URI", "sqlite:///:memory:")
    app.config.setdefault("SQLALCHEMY_BINDS", None)
    app.config.setdefault("SQLALCHEMY_TRACK_MODIFICATIONS", False)

    # Load configuration
    preset_config = load_config()

    # Config dependent section:
    if preset_config:
        try:
            _set_config_from_preset(app, preset_config)
        except Exception as e:
            print(f"Config error: {str(e)}")

    return app


def _set_config_from_preset(app: Flask, preset_config: dict):
    # Load flask info
    app.config.update(preset_config["database"])

    # Add `folders` configuration
    app.config["folders"] = preset_config.get("folders", {})

    # Ensure PubMLST configuration is included
    app.config["pubmlst"] = preset_config.get("pubmlst", {"client_id": "", "client_secret": ""})

    # Add extrapaths to config
    preset_config["folders"]["expec"] = os.path.abspath(
        os.path.join(Path(__file__).parent.parent, "unique_references/ExPEC.fsa")
    )
    # Check if release install exists
    for entry in os.listdir(get_python_lib()):
        if "microSALT-" in entry:
            preset_config["folders"]["expec"] = os.path.abspath(
                os.path.join(os.path.expandvars("$CONDA_PREFIX"), "expec/ExPEC.fsa")
            )
            break
    preset_config["folders"]["adapters"] = os.path.abspath(
        os.path.join(
            os.path.expandvars("$CONDA_PREFIX"),
            "share/trimmomatic/adapters/",
        )
    )

    # Initialize logger
    logger = initialize_logger(preset_config)

    # Create paths mentioned in config
    create_paths(preset_config, logger)

    # Integrity check database
    db_file = re.search(
        "sqlite:///(.+)",
        preset_config["database"]["SQLALCHEMY_DATABASE_URI"],
    )[1]
    check_database_integrity(db_file, logger)
