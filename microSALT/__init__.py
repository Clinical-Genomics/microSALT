import collections
import logging
import json
import os
import pathlib
import re
import subprocess
import sys

from logging import Logger
from flask import Flask
from distutils.sysconfig import get_python_lib

__version__ = "4.2.0"

app = Flask(__name__, template_folder="server/templates")
app.config.setdefault("SQLALCHEMY_DATABASE_URI", "sqlite:///:memory:")
app.config.setdefault("SQLALCHEMY_BINDS", None)
app.config.setdefault("SQLALCHEMY_TRACK_MODIFICATIONS", False)

# Keep track of microSALT installation
wd = os.path.dirname(os.path.realpath(__file__))

# Load configuration
preset_config = ""

logger: Logger | None = None

logging_levels = {
    "DEBUG": logging.DEBUG,
    "INFO": logging.INFO,
    "WARNING": logging.WARNING,
    "ERROR": logging.ERROR,
    "CRITICAL": logging.CRITICAL,
}


def setup_logger(logging_level: str) -> None:
    global logger
    logger = logging.getLogger("main_logger")
    logger.setLevel(logging_levels[logging_level])
    ch = logging.StreamHandler()
    ch.setLevel(logging_levels[logging_level])
    formatter = logging.Formatter("%(asctime)s\t%(levelname)s\t%(message)s", "%Y-%m %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)


default = os.path.join(os.environ["HOME"], ".microSALT/config.json")

if "MICROSALT_CONFIG" in os.environ:
    try:
        envvar = os.environ["MICROSALT_CONFIG"]
        with open(envvar, "r") as conf:
            preset_config = json.load(conf)
    except Exception as e:
        print("Config error: {}".format(str(e)))
        pass
elif os.path.exists(default):
    try:
        with open(os.path.abspath(default), "r") as conf:
            preset_config = json.load(conf)
    except Exception as e:
        print("Config error: {}".format(str(e)))
        pass

# Config dependent section:
if preset_config != "":
    try:
        # Load flask info
        app.config.update(preset_config["database"])

        # Add `folders` configuration
        app.config["folders"] = preset_config.get("folders", {})

        # Ensure PubMLST configuration is included

        app.config["pubmlst"] = preset_config.get("pubmlst", {"client_id": "", "client_secret": ""})

        app.config["pasteur"] = preset_config.get("pasteur", {"client_id": "", "client_secret": ""})

        # Add extrapaths to config
        preset_config["folders"]["expec"] = os.path.abspath(
            os.path.join(pathlib.Path(__file__).parent.parent, "unique_references/ExPEC.fsa")
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
        setup_logger(logging_level=logging_levels["INFO"])

        # Create paths mentioned in config
        db_file = re.search(
            "sqlite:///(.+)",
            preset_config["database"]["SQLALCHEMY_DATABASE_URI"],
        ).group(1)
        for entry in preset_config.keys():
            if entry != "_comment":
                if (
                    isinstance(preset_config[entry], str)
                    and "/" in preset_config[entry]
                    and entry not in ["genologics"]
                ):
                    if not preset_config[entry].startswith("/"):
                        sys.exit(-1)
                    unmade_fldr = os.path.abspath(preset_config[entry])
                    if not pathlib.Path(unmade_fldr).exists():
                        os.makedirs(unmade_fldr)
                        logger.info("Created path {}".format(unmade_fldr))

                # level two
                elif isinstance(preset_config[entry], collections.Mapping):
                    for thing in preset_config[entry].keys():
                        if (
                            isinstance(preset_config[entry][thing], str)
                            and "/" in preset_config[entry][thing]
                            and entry not in ["genologics"]
                        ):
                            # Special string, mangling
                            if thing == "log_file":
                                unmade_fldr = os.path.dirname(preset_config[entry][thing])
                                bash_cmd = "touch {}".format(preset_config[entry][thing])
                                proc = subprocess.Popen(bash_cmd.split(), stdout=subprocess.PIPE)
                                output, error = proc.communicate()
                            elif thing == "SQLALCHEMY_DATABASE_URI":
                                unmade_fldr = os.path.dirname(db_file)
                                bash_cmd = "touch {}".format(db_file)
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
                                logger.info("Created path {}".format(unmade_fldr))

        fh = logging.FileHandler(os.path.expanduser(preset_config["folders"]["log_file"]))
        fh.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logger.addHandler(fh)

        # Integrity check database
        cmd = "sqlite3 {0}".format(db_file)
        cmd = cmd.split()
        cmd.append("pragma integrity_check;")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        output, error = proc.communicate()
        if not "ok" in str(output):
            logger.error("Database integrity failed! Lock-state detected!")
            sys.exit(-1)

    except Exception as e:
        print("Config error: {}".format(str(e)))
        pass
