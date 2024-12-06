import base64
import hashlib
import hmac
import json
import os
import time
from datetime import datetime, timedelta
from urllib.parse import quote_plus, urlencode

from dateutil import parser
from rauth import OAuth1Session

import microSALT.utils.pubmlst.credentials as credentials

BASE_API = "https://rest.pubmlst.org"
SESSION_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "session_credentials.json")
SESSION_EXPIRATION_BUFFER = 60  # Seconds before expiration to renew


def save_session_token(token, secret, expiration_date):
    """Save session token, secret, and expiration to a JSON file."""
    session_data = {
        "token": token,
        "secret": secret,
        "expiration": expiration_date.isoformat(),
    }
    with open(SESSION_FILE, "w") as f:
        json.dump(session_data, f)
    print(f"Session token saved to {SESSION_FILE}.")


def load_session_token():
    """Load session token from file if it exists and is valid."""
    if os.path.exists(SESSION_FILE):
        with open(SESSION_FILE, "r") as f:
            session_data = json.load(f)
            expiration = parser.parse(session_data["expiration"])
            if datetime.now() < expiration - timedelta(seconds=SESSION_EXPIRATION_BUFFER):
                print("Using existing session token.")
                return session_data["token"], session_data["secret"]
    return None, None


def generate_oauth_header(url, token, token_secret):
    """Generate the OAuth1 Authorization header."""
    oauth_timestamp = str(int(time.time()))
    oauth_nonce = base64.urlsafe_b64encode(os.urandom(32)).decode("utf-8").strip("=")
    oauth_signature_method = "HMAC-SHA1"
    oauth_version = "1.0"

    oauth_params = {
        "oauth_consumer_key": credentials.CLIENT_ID,
        "oauth_token": token,
        "oauth_signature_method": oauth_signature_method,
        "oauth_timestamp": oauth_timestamp,
        "oauth_nonce": oauth_nonce,
        "oauth_version": oauth_version,
    }

    # Create the signature base string
    params_encoded = urlencode(sorted(oauth_params.items()))
    base_string = f"GET&{quote_plus(url)}&{quote_plus(params_encoded)}"
    signing_key = f"{credentials.CLIENT_SECRET}&{token_secret}"

    # Sign the base string
    hashed = hmac.new(signing_key.encode("utf-8"), base_string.encode("utf-8"), hashlib.sha1)
    oauth_signature = base64.b64encode(hashed.digest()).decode("utf-8")

    # Add the signature
    oauth_params["oauth_signature"] = oauth_signature

    # Construct the Authorization header
    auth_header = "OAuth " + ", ".join(
        [f'{quote_plus(k)}="{quote_plus(v)}"' for k, v in oauth_params.items()]
    )
    return auth_header


def get_new_session_token():
    """Request a new session token using client credentials."""
    print("Fetching a new session token...")
    db = "pubmlst_neisseria_seqdef"
    url = f"{BASE_API}/db/{db}/oauth/get_session_token"

    session = OAuth1Session(
        consumer_key=credentials.CLIENT_ID,
        consumer_secret=credentials.CLIENT_SECRET,
        access_token=credentials.ACCESS_TOKEN,
        access_token_secret=credentials.ACCESS_SECRET,
    )

    response = session.get(url, headers={"User-Agent": "BIGSdb downloader"})
    if response.status_code == 200:
        token_data = response.json()
        session_token = token_data["oauth_token"]
        session_secret = token_data["oauth_token_secret"]
        expiration_time = datetime.now() + timedelta(hours=12)  # 12-hour validity
        save_session_token(session_token, session_secret, expiration_time)
        return session_token, session_secret
    else:
        print(f"Error: {response.status_code} - {response.text}")
        return None, None
