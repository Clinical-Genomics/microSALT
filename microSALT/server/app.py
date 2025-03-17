
import logging
import os
from pathlib import Path
import re
import subprocess
import sys

from flask import Flask
from microsalt.server.views import bp
from microsalt.store.orm_models import db
logger = logging.getLogger("main_logger")

APP: Flask = None

def initialize_app(config: dict) -> None:
    global APP
    APP = Flask(__name__, template_folder="templates")
    APP.config.setdefault("SQLALCHEMY_DATABASE_URI", "sqlite:///:memory:")
    APP.config.setdefault("SQLALCHEMY_BINDS", {})
    APP.config.setdefault("SQLALCHEMY_TRACK_MODIFICATIONS", False)
    _setup_app_from_config(APP, config)
    db.init_app(APP)
    APP.register_blueprint(bp)

def get_app() -> Flask:
    if not APP:
        raise ValueError("App not initialized")
    return APP

def create_path(path: str, logger: logging.Logger):
    if not Path(path).exists():
        os.makedirs(path)
        logger.info(f"Created path {path}")
    else:
        logger.info(f"Path {path} already exists")


def handle_config_entry(entry, value, logger):
    
    if isinstance(value, str) and "/" in value and entry not in ["genologics"]:
        if os.path.abspath(value) == "/path":
            logger.error(f"Path for {entry} is not set.")
        if not value.startswith("/"):
            sys.exit(-1)
        create_path(os.path.abspath(value), logger)
    elif isinstance(value, dict):
        for sub_entry, sub_value in value.items():
            handle_config_entry(sub_entry, sub_value, logger)


def create_paths(config, logger):

    for entry, value in config.items():
        if entry == "_comment":
            continue
        handle_config_entry(entry, value, logger)


def check_database_integrity(db_file, logger):
    if not os.path.exists(db_file):
        logger.error(f"Database file {db_file} does not exist.")
        sys.exit(-1)
    cmd = ["sqlite3", db_file, "PRAGMA integrity_check;"]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    if process.returncode != 0:
        logger.error(f"Database integrity check failed with return code {process.returncode}.")
        logger.error(f"Error: {error.decode('utf-8')}")
        sys.exit(-1)

    if "ok" not in output.decode("utf-8"):
        logger.error("Database integrity failed! Lock-state detected!")
        sys.exit(-1)

    logger.info("Database integrity check passed.")

def _setup_app_from_config(app: Flask, config: dict):
    app.config.update(config["database"])
    app.config["folders"] = config.get("folders", {})
    app.config["pubmlst"] = config.get("pubmlst", {"client_id": "", "client_secret": ""})
    app.config["threshold"] = config.get("threshold", {})
    app.config["verified_organisms_regex"] = config["regex"]["verified_organisms"]
    app.config["mlst_novel_id_threshold"] = config["threshold"]["mlst_novel_id"]
    app.config["mlst_id_threshold"] = config["threshold"]["mlst_id"]
    app.config["mlst_span_threshold"] = config["threshold"]["mlst_span"]
    app.config["motif_id_threshold"] = config["threshold"]["motif_id"]
    app.config["motif_span_threshold"] = config["threshold"]["motif_span"]
    
    
    
    # Create paths mentioned in config
    create_paths(app.config["folders"], logger)
    
    # Integrity check database
    db_file = re.search(
        "sqlite:///(.+)",
        config["database"]["SQLALCHEMY_DATABASE_URI"],
    )[1]
    check_database_integrity(db_file, logger)

if __name__ == "__main__":
    
    from microSALT.utils.config_loader import load_config
    config = load_config('/Users/karlnyren/github/microSALT/tests/testdata/config.json')
    initialize_app(config)
    app = get_app()
    app.run(host="127.0.0.1", port=5001, debug=True, use_reloader=False)
