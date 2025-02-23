import json
import os
from datetime import datetime, timedelta
from dateutil import parser
from rauth import OAuth1Session
from microSALT import logger
from microSALT.utils.pubmlst.helpers import BASE_API, save_session_token, load_auth_credentials, get_path, folders_config, credentials_path_key, pubmlst_session_credentials_file_name 
from microSALT.utils.pubmlst.exceptions import (
    PUBMLSTError,
    SessionTokenRequestError,
    SessionTokenResponseError,
)

session_token_validity = 12  # 12-hour validity
session_expiration_buffer = 60  # 60-second buffer

def get_new_session_token(db: str):
    """Request a new session token using all credentials for a specific database."""
    logger.debug("Fetching a new session token for database '{db}'...")

    try:
        consumer_key, consumer_secret, access_token, access_secret = load_auth_credentials()

        url = f"{BASE_API}/db/{db}/oauth/get_session_token"

        session = OAuth1Session(
            consumer_key=consumer_key,
            consumer_secret=consumer_secret,
            access_token=access_token,
            access_token_secret=access_secret,
        )

        response = session.get(url, headers={"User-Agent": "BIGSdb downloader"})
        logger.debug("Response Status Code: {status_code}")

        if response.ok:
            try:
                token_data = response.json()
                session_token = token_data.get("oauth_token")
                session_secret = token_data.get("oauth_token_secret")

                if not session_token or not session_secret:
                    raise SessionTokenResponseError(
                        db, "Missing 'oauth_token' or 'oauth_token_secret' in response."
                    )

                expiration_time = datetime.now() + timedelta(hours=session_token_validity)

                save_session_token(db, session_token, session_secret, expiration_time)
                return session_token, session_secret

            except (ValueError, KeyError) as e:
                raise SessionTokenResponseError(db, f"Invalid response format: {str(e)}")
        else:
            raise SessionTokenRequestError(
                db, response.status_code, response.text
            )

    except PUBMLSTError as e:
        logger.error(f"Error during token fetching: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise PUBMLSTError(f"Unexpected error while fetching session token for database '{db}': {e}")

def load_session_credentials(db: str):
    """Load session token from file for a specific database."""
    try:
        credentials_file = os.path.join(
            get_path(folders_config, credentials_path_key),
            pubmlst_session_credentials_file_name
        )

        if not os.path.exists(credentials_file):
            logger.debug("Session file does not exist. Fetching a new session token.")
            return get_new_session_token(db)

        with open(credentials_file, "r") as f:
            try:
                all_sessions = json.load(f)
            except json.JSONDecodeError as e:
                raise SessionTokenResponseError(db, f"Failed to parse session file: {str(e)}")

        db_session_data = all_sessions.get("databases", {}).get(db)
        if not db_session_data:
            logger.debug(f"No session token found for database '{db}'. Fetching a new session token.")
            return get_new_session_token(db)

        expiration = parser.parse(db_session_data.get("expiration", ""))
        if datetime.now() < expiration - timedelta(seconds=session_expiration_buffer):
            logger.debug(f"Using existing session token for database '{db}'.")
            session_token = db_session_data.get("token")
            session_secret = db_session_data.get("secret")

            return session_token, session_secret

        logger.debug(f"Session token for database '{db}' has expired. Fetching a new session token.")
        return get_new_session_token(db)

    except PUBMLSTError as e:
        logger.error(f"PUBMLST-specific error occurred: {e}")
        raise
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise PUBMLSTError(f"Unexpected error while loading session token for database '{db}': {e}")

