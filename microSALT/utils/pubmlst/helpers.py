import os
import requests
import base64
import hashlib
import hmac
import json
import time
from pathlib import Path
from urllib.parse import quote_plus, urlencode

from werkzeug.exceptions import NotFound

from microSALT import app, logger
from microSALT.utils.pubmlst.constants import Encoding, url_map
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
pubmlst_config = app.config["pubmlst"]
folders_config = app.config["folders"]


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

def get_db_type_capabilities(db_name: str) -> dict:
    """
    Determine whether the database is of type 'isolate' or 'sequence definition (seqdef)'.
    This is inferred by inspecting metadata for known capabilities.
    """
    url = f"{BASE_API}/db/{db_name}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        metadata = response.json()
        capabilities = {
            "has_isolates": metadata.get("has_isolates", False),
            "has_projects": metadata.get("has_projects", False),
            "has_fields": metadata.get("has_fields", False),
        }
        return capabilities
    except Exception as e:
        raise PUBMLSTError(f"Failed to get DB type capabilities for {db_name}: {e}") from e


def should_skip_endpoint(endpoint: str, capabilities: dict) -> bool:
    """
    Return True if the endpoint call should be skipped due to being incompatible
    with the database's declared capabilities.
    """
    # Handle isolate-related endpoints
    if "/isolates" in endpoint and not capabilities.get("has_isolates", False):
        return True
    if "/projects" in endpoint and not capabilities.get("has_projects", False):
        return True
    if "/fields" in endpoint and not capabilities.get("has_fields", False):
        return True
    return False
