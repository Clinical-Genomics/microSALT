import json
import os
from datetime import datetime, timedelta
from pathlib import Path
from dateutil import parser
from rauth import OAuth1Session
from microSALT import app, logger
from microSALT.utils.pubmlst.helpers import get_credentials_file_path, BASE_API, load_credentials, generate_oauth_header

SESSION_EXPIRATION_BUFFER = 60  # Seconds before expiration to renew

pubmlst_config = app.config["pubmlst"]
credentials_files_path = get_credentials_file_path(pubmlst_config)

# Ensure the directory exists
credentials_files_path.mkdir(parents=True, exist_ok=True)

CREDENTIALS_FILE = os.path.join(credentials_files_path, "PUBMLST_credentials.py")
SESSION_FILE = os.path.join(credentials_files_path, "PUBMLST_session_credentials.json")


def save_session_token(db, token, secret, expiration_date):
    """Save session token, secret, and expiration to a JSON file for the specified database."""
    session_data = {
        "token": token,
        "secret": secret,
        "expiration": expiration_date.isoformat(),
    }

    # Load existing sessions if available
    if os.path.exists(SESSION_FILE):
        with open(SESSION_FILE, "r") as f:
            all_sessions = json.load(f)
    else:
        all_sessions = {}

    # Ensure 'databases' key exists
    if "databases" not in all_sessions:
        all_sessions["databases"] = {}

    # Update the session token for the specific database
    all_sessions["databases"][db] = session_data

    # Save back to file
    with open(SESSION_FILE, "w") as f:
        json.dump(all_sessions, f, indent=4)
    logger.info(f"Session token for '{db}' saved to {SESSION_FILE}.")


def load_session_token(db):
    """Load session token from file for a specific database if it exists and is valid."""
    if not os.path.exists(SESSION_FILE):
        logger.info("Session file does not exist.")
        return None, None

    with open(SESSION_FILE, "r") as f:
        all_sessions = json.load(f)

    # Check if the database entry exists
    db_session_data = all_sessions.get("databases", {}).get(db)
    if not db_session_data:
        logger.info(f"No session token found for database '{db}'.")
        return None, None

    expiration = parser.parse(db_session_data["expiration"])
    if datetime.now() < expiration - timedelta(seconds=SESSION_EXPIRATION_BUFFER):
        logger.info(f"Using existing session token for database '{db}'.")
        return db_session_data["token"], db_session_data["secret"]
    else:
        logger.info(f"Session token for database '{db}' has expired.")
        return None, None


def get_new_session_token(db="pubmlst_test_seqdef"):
    """Request a new session token using all credentials for a specific database."""
    logger.info(f"Fetching a new session token for database '{db}'...")
    client_id, client_secret, access_token, access_secret = load_credentials()
    url = f"{BASE_API}/db/{db}/oauth/get_session_token"

    # Create an OAuth1Session with all credentials
    session = OAuth1Session(
        consumer_key=client_id,
        consumer_secret=client_secret,
        access_token=access_token,
        access_token_secret=access_secret,
    )

    try:
        response = session.get(url, headers={"User-Agent": "BIGSdb downloader"})
        logger.info(f"Response Status Code: {response.status_code}")


        if response.status_code == 200:
            token_data = response.json()
            session_token = token_data["oauth_token"]
            session_secret = token_data["oauth_token_secret"]
            expiration_time = datetime.now() + timedelta(hours=12)  # 12-hour validity
            save_session_token(db, session_token, session_secret, expiration_time)
            return session_token, session_secret
        else:
            raise ValueError(f"Error fetching session token: {response.status_code} - {response.text}")
    except Exception as e:
        logger.error(f"Error during token fetching: {e}")
        raise
