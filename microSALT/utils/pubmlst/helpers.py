import os
import requests
import base64
import hashlib
import json
import hmac
import time
from pathlib import Path
from urllib.parse import quote_plus, urlencode
from microSALT import app, logger
from microSALT.utils.pubmlst.exceptions import PathResolutionError, CredentialsFileNotFound, InvalidCredentials, PubMLSTError, SaveSessionError, InvalidURLError
from microSALT.utils.pubmlst.constants import CREDENTIALS_KEY, Encoding, URL_MAPS

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


def get_service_config(service: str):
    """
    Get the configuration for the specified service (e.g., 'pubmlst' or 'pasteur').

    :param service: The name of the service ('pubmlst' or 'pasteur').
    :return: A dictionary containing the configuration for the service.
    """
    services = {
        "pubmlst": {
            "base_web": "https://pubmlst.org/bigsdb",
            "base_api": "https://rest.pubmlst.org",
            "base_api_host": "rest.pubmlst.org",
            "database": "pubmlst_test_seqdef",
            "auth_credentials_file_name": "pubmlst_credentials.env",
            "session_credentials_file_name": "pubmlst_session_credentials.json",
            "config": app.config["pubmlst"],
        },
        "pasteur": {
            "base_web": "https://bigsdb.pasteur.fr/cgi-bin/bigsdb/bigsdb.pl",
            "base_api": "https://bigsdb.pasteur.fr/api",
            "base_api_host": "bigsdb.pasteur.fr",
            "auth_credentials_file_name": "pasteur_credentials.env",
            "session_credentials_file_name": "pasteur_session_credentials.json",
            "config": app.config["pasteur"],
        },
    }

    if service not in services:
        raise ValueError(f"Unknown service: {service}")

    return services[service]


def get_service_by_url(url: str):
    """
    Get the client name based on the provided URL.

    :param url: The URL to check.
    :return: The name of the client ('pubmlst' or 'pasteur').
    """
    for service in ["pubmlst", "pasteur"]:
        if url.startswith(get_service_config(service)["base_api"]):
            return service
    raise ValueError(f"Unknown client for URL: {url}")


def get_url_map(service: str):
    """
    Get the URL map for the specified service.

    :param service: The name of the service ('pubmlst' or 'pasteur').
    :return: The URL map for the service.
    """
    url_map = URL_MAPS.get(service, None)
    if url_map is None:
        raise ValueError(f"Unknown service: {service}")
    return url_map


def load_auth_credentials(service: str):
    """
    Load client ID, client secret, access token, and access secret from the credentials file for the specified service.

    :param service: The name of the service ('pubmlst' or 'pasteur').
    :return: A tuple containing the credentials (consumer_key, consumer_secret, access_token, access_secret).
    """
    try:
        service_config = get_service_config(service)
        credentials_file = os.path.join(
            get_path(folders_config, CREDENTIALS_KEY),
            service_config["auth_credentials_file_name"],
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
    except PubMLSTError as e:
        logger.error(f"Unexpected error in load_{service}_credentials: {e}")
        raise
    except Exception as e:
        raise PubMLSTError(f"An unexpected error occurred while loading {service} credentials: {e}")
    

def generate_oauth_header(url: str, oauth_consumer_key: str, oauth_consumer_secret: str, oauth_token: str, oauth_token_secret: str):
    """Generate the OAuth1 Authorization header."""
    oauth_timestamp = str(int(time.time()))
    oauth_nonce = base64.urlsafe_b64encode(os.urandom(32)).decode(Encoding.UTF8.value).strip("=")
    oauth_signature_method = "HMAC-SHA1"
    oauth_version = "1.0"

    oauth_params = {
        "oauth_consumer_key": oauth_consumer_key,
        "oauth_token": oauth_token,
        "oauth_signature_method": oauth_signature_method,
        "oauth_timestamp": oauth_timestamp,
        "oauth_nonce": oauth_nonce,
        "oauth_version": oauth_version,
    }    

    params_encoded = urlencode(sorted(oauth_params.items()))
    base_string = f"GET&{quote_plus(url)}&{quote_plus(params_encoded)}"
    signing_key = f"{oauth_consumer_secret}&{oauth_token_secret}"

    hashed = hmac.new(signing_key.encode(Encoding.UTF8.value), base_string.encode(Encoding.UTF8.value), hashlib.sha1)
    oauth_signature = base64.b64encode(hashed.digest()).decode(Encoding.UTF8.value)

    oauth_params["oauth_signature"] = oauth_signature

    auth_header = "OAuth " + ", ".join(
        [f'{quote_plus(k)}="{quote_plus(v)}"' for k, v in oauth_params.items()]
    )
    return auth_header


def save_session_token(service: str, db: str, token: str, secret: str, expiration_date: str):
    """
    Save session token, secret, and expiration to a JSON file for the specified service and database.

    :param service: The name of the service ('pubmlst' or 'pasteur').
    :param db: The database name.
    :param token: The session token.
    :param secret: The session secret.
    :param expiration_date: The expiration date of the session token.
    """
    try:
        service_config = get_service_config(service)
        session_file = os.path.join(
            get_path(folders_config, CREDENTIALS_KEY),
            service_config["session_credentials_file_name"],
        )

        session_data = {
            "token": token,
            "secret": secret,
            "expiration": expiration_date.isoformat(),
        }

        if os.path.exists(session_file):
            with open(session_file, "r") as f:
                all_sessions = json.load(f)
        else:
            all_sessions = {}

        if "databases" not in all_sessions:
            all_sessions["databases"] = {}

        all_sessions["databases"][db] = session_data

        with open(session_file, "w") as f:
            json.dump(all_sessions, f, indent=4)

        logger.debug(f"Session token for {service} database '{db}' saved to '{session_file}'.")

    except (IOError, OSError) as e:
        raise SaveSessionError(db, f"I/O error: {e}")
    except ValueError as e:
        raise SaveSessionError(db, f"Invalid data format: {e}")
    except Exception as e:
        raise SaveSessionError(db, f"Unexpected error: {e}")


def get_db_type_capabilities(base_api: str, db_name: str) -> dict:
    """
    Determine whether the database is of type 'isolate' or 'sequence definition (seqdef)'.
    This is inferred by inspecting metadata for known capabilities.
    """
    url = f"{base_api}/db/{db_name}"
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
        raise PubMLSTError(f"Failed to get DB type capabilities for {db_name}: {e}") from e


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
