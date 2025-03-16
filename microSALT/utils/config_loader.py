import json
import os
import importlib.util
from pathlib import Path


def load_config(config_path: str = None) -> dict:
    config = {}

    if config_path:
        with open(config_path, "r") as conf:
            config = json.load(conf)

    # Get the path to the microSALT package
    microsalt_spec = importlib.util.find_spec("microSALT")
    if microsalt_spec and microsalt_spec.origin:
        microsalt_dir_path = Path(microsalt_spec.origin).parent
        config["folders"]["expec"] = Path(microsalt_dir_path, "unique_references/ExPEC.fsa")

    return config
