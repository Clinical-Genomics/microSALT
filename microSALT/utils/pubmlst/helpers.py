import logging
import json
import os

from flask import Flask
from pathlib import Path

from werkzeug.exceptions import NotFound

from microSALT.server.app import get_app
from microSALT.utils.pubmlst.constants import url_map
from microSALT.utils.pubmlst.exceptions import (
    CredentialsFileNotFound,
    InvalidCredentials,
    InvalidURLError,
    PathResolutionError,
    PUBMLSTError,
    SaveSessionError,
)

BASE_WEB = "https://pubmlst.org/bigsdb"
BASE_API = "https://rest.pubmlst.org"
BASE_API_HOST = "rest.pubmlst.org"

credentials_path_key = "pubmlst_credentials"
pubmlst_auth_credentials_file_name = "pubmlst_credentials.env"
pubmlst_session_credentials_file_name = "pubmlst_session_credentials.json"

logger = logging.getLogger(__name__)


def get_folders_config():
    """Get the folders configuration from the application."""
    app: Flask = get_app()
    return app.config["folders"]


def get_path(config, config_key: str):
    """Get and expand the file path from the configuration."""
    try:
        path = config.get(config_key)
        if not path:
            raise PathResolutionError(config_key)

        path = os.path.expandvars(path)
        path = os.path.expanduser(path)

        return Path(path).resolve()

    except Exception as e:
        raise PathResolutionError(config_key) from e


def load_auth_credentials():
    """Load client ID, client secret, access token, and access secret from credentials file."""
    folders_config = get_folders_config()
    try:
        credentials_file = os.path.join(
            get_path(folders_config, credentials_path_key), pubmlst_auth_credentials_file_name
        )

        if not os.path.exists(credentials_file):
            raise CredentialsFileNotFound(credentials_file)

        credentials = {}
        with open(credentials_file, "r") as f:
            exec(f.read(), credentials)

        consumer_key = credentials.get("CLIENT_ID", "").strip()
        consumer_secret = credentials.get("CLIENT_SECRET", "").strip()
        access_token = credentials.get("ACCESS_TOKEN", "").strip()
        access_secret = credentials.get("ACCESS_SECRET", "").strip()

        missing_fields = []
        if not consumer_key:
            missing_fields.append("CLIENT_ID")
        if not consumer_secret:
            missing_fields.append("CLIENT_SECRET")
        if not access_token:
            missing_fields.append("ACCESS_TOKEN")
        if not access_secret:
            missing_fields.append("ACCESS_SECRET")

        if missing_fields:
            raise InvalidCredentials(missing_fields)

        return consumer_key, consumer_secret, access_token, access_secret

    except CredentialsFileNotFound:
        raise
    except InvalidCredentials:
        raise
    except PUBMLSTError as e:
        logger.error(f"Unexpected error in load_credentials: {e}")
        raise
    except Exception as e:
        raise PUBMLSTError("An unexpected error occurred while loading credentials: {e}")


def save_session_token(db: str, token: str, secret: str, expiration_date: str):
    """Save session token, secret, and expiration to a JSON file for the specified database."""
    folders_config = get_folders_config()
    try:
        session_data = {
            "token": token,
            "secret": secret,
            "expiration": expiration_date.isoformat(),
        }

        credentials_file = os.path.join(
            get_path(folders_config, credentials_path_key), pubmlst_session_credentials_file_name
        )

        if os.path.exists(credentials_file):
            with open(credentials_file, "r") as f:
                all_sessions = json.load(f)
        else:
            all_sessions = {}

        if "databases" not in all_sessions:
            all_sessions["databases"] = {}

        all_sessions["databases"][db] = session_data

        with open(credentials_file, "w") as f:
            json.dump(all_sessions, f, indent=4)

        logger.debug(f"Session token for database '{db}' saved to '{credentials_file}'.")
    except (IOError, OSError) as e:
        raise SaveSessionError(db, f"I/O error: {e}")
    except ValueError as e:
        raise SaveSessionError(db, f"Invalid data format: {e}")
    except Exception as e:
        raise SaveSessionError(db, f"Unexpected error: {e}")


def parse_pubmlst_url(url: str):
    """
    Match a URL against the URL map and return extracted parameters.
    """
    adapter = url_map.bind("")
    parsed_url = url.split(BASE_API_HOST)[-1]
    try:
        endpoint, values = adapter.match(parsed_url)
        return {"endpoint": endpoint, **values}
    except NotFound:
        raise InvalidURLError(url)
