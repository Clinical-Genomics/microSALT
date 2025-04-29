import json
import os
from pathlib import Path

from werkzeug.routing import Map, Rule

from microSALT import app, logger

from microSALT.utils.pubmlst.exceptions import (
    CredentialsFileNotFound,
    InvalidCredentials,
    PathResolutionError,
    PUBMLSTError,
    SaveSessionError,
)
from microSALT.utils.pubmlst.constants import CREDENTIALS_KEY, URL_MAPS

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
    except PUBMLSTError as e:
        logger.error(f"Unexpected error in load_{service}_credentials: {e}")
        raise
    except Exception as e:
        raise PUBMLSTError(f"An unexpected error occurred while loading {service} credentials: {e}")


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
