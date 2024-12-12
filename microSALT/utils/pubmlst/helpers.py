import os
import json
import base64
import hashlib
import hmac
import time
from pathlib import Path
from urllib.parse import quote_plus, urlencode
import requests
from datetime import datetime, timedelta
from dateutil import parser
from microSALT import app, logger

BASE_WEB = {
    "PubMLST": "https://pubmlst.org/bigsdb",
}

BASE_API_DICT = {
    "PubMLST": "https://rest.pubmlst.org",
}

BASE_API = "https://rest.pubmlst.org"  # Used by authentication and other modules

def get_credentials_file_path(pubmlst_config):
    """Get and expand the credentials file path from the configuration."""
    # Retrieve the path from config or use current working directory if not set
    path = pubmlst_config.get("credentials_files_path", os.getcwd())
    # Expand environment variables like $HOME
    path = os.path.expandvars(path)
    # Expand user shortcuts like ~
    path = os.path.expanduser(path)
    return Path(path).resolve()

def load_credentials():
    """Load client ID, client secret, access token, and access secret from credentials file."""
    pubmlst_config = app.config["pubmlst"]
    credentials_files_path = get_credentials_file_path(pubmlst_config)
    credentials_file = os.path.join(credentials_files_path, "PUBMLST_credentials.py")

    if not os.path.exists(credentials_file):
        raise FileNotFoundError(
            f"Credentials file not found: {credentials_file}. "
            "Please generate it using get_credentials.py."
        )
    credentials = {}
    with open(credentials_file, "r") as f:
        exec(f.read(), credentials)

    client_id = credentials.get("CLIENT_ID", "").strip()
    client_secret = credentials.get("CLIENT_SECRET", "").strip()
    access_token = credentials.get("ACCESS_TOKEN", "").strip()
    access_secret = credentials.get("ACCESS_SECRET", "").strip()

    if not (client_id and client_secret and access_token and access_secret):
        raise ValueError(
            "Invalid credentials: All fields (CLIENT_ID, CLIENT_SECRET, ACCESS_TOKEN, ACCESS_SECRET) must be non-empty. "
            "Please regenerate the credentials file using get_credentials.py."
        )
    return client_id, client_secret, access_token, access_secret

def generate_oauth_header(url, token, token_secret):
    """Generate the OAuth1 Authorization header."""
    client_id, client_secret, _, _ = load_credentials()
    oauth_timestamp = str(int(time.time()))
    oauth_nonce = base64.urlsafe_b64encode(os.urandom(32)).decode("utf-8").strip("=")
    oauth_signature_method = "HMAC-SHA1"
    oauth_version = "1.0"

    oauth_params = {
        "oauth_consumer_key": client_id,
        "oauth_token": token,
        "oauth_signature_method": oauth_signature_method,
        "oauth_timestamp": oauth_timestamp,
        "oauth_nonce": oauth_nonce,
        "oauth_version": oauth_version,
    }

    params_encoded = urlencode(sorted(oauth_params.items()))
    base_string = f"GET&{quote_plus(url)}&{quote_plus(params_encoded)}"
    signing_key = f"{client_secret}&{token_secret}"

    hashed = hmac.new(signing_key.encode("utf-8"), base_string.encode("utf-8"), hashlib.sha1)
    oauth_signature = base64.b64encode(hashed.digest()).decode("utf-8")

    oauth_params["oauth_signature"] = oauth_signature

    auth_header = "OAuth " + ", ".join(
        [f'{quote_plus(k)}="{quote_plus(v)}"' for k, v in oauth_params.items()]
    )
    return auth_header

def validate_session_token(session_token, session_secret):
    """Ensure session token and secret are valid."""
    if not session_token or not session_secret:
        raise ValueError("Session token or secret is missing. Please authenticate first.")

def fetch_paginated_data(url, session_token, session_secret):
    """Fetch paginated data using the session token and secret."""
    validate_session_token(session_token, session_secret)

    results = []
    while url:
        headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
        response = requests.get(url, headers=headers)

        logger.debug(f"Fetching URL: {url}")
        logger.debug(f"Response Status Code: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            results.extend(data.get("profiles", []))
            url = data.get("paging", {}).get("next", None)  # Get the next page URL if available
        else:
            raise ValueError(
                f"Failed to fetch data. URL: {url}, Status Code: {response.status_code}, "
                f"Response: {response.text}"
            )
    return results
