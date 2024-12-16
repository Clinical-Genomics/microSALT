import os
import base64
import hashlib
import json
import hmac
import time
from pathlib import Path
from urllib.parse import quote_plus, urlencode
from microSALT import app, logger
from microSALT.utils.pubmlst.exceptions import PUBMLSTError, PathResolutionError, CredentialsFileNotFound, InvalidCredentials, SaveSessionError
from microSALT.utils.pubmlst.constants import Encoding

BASE_WEB = "https://pubmlst.org/bigsdb"
BASE_API = "https://rest.pubmlst.org" 
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
            get_path(folders_config, credentials_path_key),
            pubmlst_auth_credentials_file_name
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

def save_session_token(db: str, token: str, secret: str, expiration_date: str):
    """Save session token, secret, and expiration to a JSON file for the specified database."""
    try:
        session_data = {
            "token": token,
            "secret": secret,
            "expiration": expiration_date.isoformat(),
        }

        credentials_file = os.path.join(
            get_path(folders_config, credentials_path_key),
            pubmlst_session_credentials_file_name
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

        logger.debug(
            f"Session token for database '{db}' saved to '{credentials_file}'."
        )
    except (IOError, OSError) as e:
        raise SaveSessionError(db, f"I/O error: {e}")
    except ValueError as e:
        raise SaveSessionError(db, f"Invalid data format: {e}")
    except Exception as e:
        raise SaveSessionError(db, f"Unexpected error: {e}")
