import collections
import logging
import json
import os
import pathlib
import re
import subprocess
import sys

from flask import Flask
from distutils.sysconfig import get_python_lib

__version__ = "4.0.0"

app = Flask(__name__, template_folder="server/templates")
app.config.setdefault("SQLALCHEMY_DATABASE_URI", "sqlite:///:memory:")
app.config.setdefault("SQLALCHEMY_BINDS", None)
app.config.setdefault("SQLALCHEMY_TRACK_MODIFICATIONS", False)

# Reusable function for resolving paths
def resolve_path(path):
    """Resolve environment variables, user shortcuts, and absolute paths."""
    if path:
        path = os.path.expandvars(path)  # Expand environment variables like $HOME
        path = os.path.expanduser(path)  # Expand user shortcuts like ~
        path = os.path.abspath(path)    # Convert to an absolute path
        return path
    return path

# Function to create directories if they do not exist
def ensure_directory(path, logger=None):
    """Ensure a directory exists; create it if missing."""
    try:
        if path and not pathlib.Path(path).exists():
            os.makedirs(path, exist_ok=True)
            if logger:
                logger.info("Created path {}".format(path))
    except Exception as e:
        if logger:
            logger.error("Failed to create path {}: {}".format(path, e))
        raise

# Initialize logger
logger = logging.getLogger("main_logger")
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))
logger.addHandler(ch)

# Keep track of microSALT installation
wd = os.path.dirname(os.path.realpath(__file__))

# Load configuration
preset_config = ""
default_config_path = resolve_path("$HOME/.microSALT/config.json")

if "MICROSALT_CONFIG" in os.environ:
    try:
        envvar = os.environ["MICROSALT_CONFIG"]
        with open(envvar, "r") as conf:
            preset_config = json.load(conf)
    except Exception as e:
        logger.error("Config error: {}".format(e))
elif os.path.exists(default_config_path):
    try:
        with open(default_config_path, "r") as conf:
            preset_config = json.load(conf)
    except Exception as e:
        logger.error("Config error: {}".format(e))

# Config dependent section:
if preset_config:
    try:
        # Load Flask info
        app.config.update(preset_config["database"])

        # Add extrapaths to config
        preset_config["folders"]["expec"] = resolve_path(
            os.path.join(pathlib.Path(__file__).parent.parent, "unique_references/ExPEC.fsa")
        )

        # Check if release install exists
        for entry in os.listdir(get_python_lib()):
            if "microSALT-" in entry:
                preset_config["folders"]["expec"] = resolve_path(
                    os.path.join(os.path.expandvars("$CONDA_PREFIX"), "expec/ExPEC.fsa")
                )
                break

        preset_config["folders"]["adapters"] = resolve_path(
            os.path.join(os.path.expandvars("$CONDA_PREFIX"), "share/trimmomatic/adapters/")
        )

        # Load pubmlst configuration
        if "pubmlst" not in preset_config:
            raise KeyError("Missing 'pubmlst' section in configuration file.")
        pubmlst_config = preset_config["pubmlst"]

        # Set and resolve credentials file path
        credentials_files_path = resolve_path(pubmlst_config.get("credentials_files_path", "$HOME/.microSALT"))
        pubmlst_config["credentials_files_path"] = credentials_files_path

        # Ensure the credentials directory exists
        ensure_directory(credentials_files_path, logger)

        # Update the app configuration
        app.config["pubmlst"] = pubmlst_config

        # Log the resolved credentials file path
        logger.info("PubMLST configuration loaded with credentials_files_path: {}".format(credentials_files_path))

    except KeyError as e:
        logger.error("Configuration error: {}".format(e))
        sys.exit(1)
    except Exception as e:
        logger.error("Unexpected error: {}".format(e))
        sys.exit(1)
