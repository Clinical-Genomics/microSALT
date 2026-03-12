import os
import sys

from requests_oauthlib import OAuth1Session

from microSALT.config import (
    MicroSALTConfig,
)
from microSALT.utils.pubmlst.constants import CREDENTIALS_KEY
from microSALT.utils.pubmlst.helpers import get_path, get_service_config


def validate_credentials(client_id, client_secret):
    """Ensure client_id and client_secret are not empty."""
    if not client_id or not client_id.strip():
        raise ValueError("Invalid CLIENT_ID: It must not be empty.")
    if not client_secret or not client_secret.strip():
        raise ValueError("Invalid CLIENT_SECRET: It must not be empty.")


def get_new_access_token(
    client_id, client_secret, db: str, base_api: str, base_web: str
) -> tuple[str, str]:
    """Obtain a new access token and secret."""
    # Step 1: fetch request token
    oauth = OAuth1Session(client_id, client_secret=client_secret, callback_uri="oob")
    response = oauth.fetch_request_token(f"{base_api}/db/{db}/oauth/get_request_token")
    if not response:
        print("Error obtaining request token.")
        sys.exit(1)
    request_token = response["oauth_token"]
    request_secret = response["oauth_token_secret"]

    print(
        "Please log in using your user account at "
        f"{base_web}?db={db}&page=authorizeClient&oauth_token={request_token} "
        "using a web browser to obtain a verification code."
    )
    verifier = input("Please enter verification code: ")

    # Step 2: exchange for access token
    oauth = OAuth1Session(
        client_id,
        client_secret=client_secret,
        resource_owner_key=request_token,
        resource_owner_secret=request_secret,
        verifier=verifier,
    )
    access_data = oauth.fetch_access_token(f"{base_api}/db/{db}/oauth/get_access_token")
    if not access_data:
        print("Error obtaining access token.")
        sys.exit(1)
    return access_data["oauth_token"], access_data["oauth_token_secret"]


def save_to_credentials_py(
    client_id, client_secret, access_token, access_secret, credentials_path, credentials_file
) -> None:
    """Save tokens in the credentials.py file."""
    credentials_path.mkdir(parents=True, exist_ok=True)

    with open(credentials_file, "w") as f:
        f.write(f'CLIENT_ID = "{client_id}"\n')
        f.write(f'CLIENT_SECRET = "{client_secret}"\n')
        f.write(f'ACCESS_TOKEN = "{access_token}"\n')
        f.write(f'ACCESS_SECRET = "{access_secret}"\n')
    print(f"Tokens saved to {credentials_file}")


def main(service: str, config: MicroSALTConfig, species: str | None = None):
    try:
        service_config = get_service_config(service, pubmlst=config.pubmlst, pasteur=config.pasteur)
        bigsd_config = service_config["config"]
        client_id = bigsd_config.client_id
        client_secret = bigsd_config.client_secret
        validate_credentials(client_id, client_secret)

        # Determine the database
        if service == "pubmlst":
            db = service_config["database"]
        elif service == "pasteur":
            if not species:
                raise ValueError("For the 'pasteur' service, you must provide a species.")
            db = f"pubmlst_{species}_seqdef"
        else:
            raise ValueError(f"Unknown service: {service}")

        credentials_path = get_path(folders=config.folders, config_key=CREDENTIALS_KEY)
        credentials_file = os.path.join(
            credentials_path, service_config.get("auth_credentials_file_name")
        )

        access_token, access_secret = get_new_access_token(
            client_id=client_id,
            client_secret=client_secret,
            db=db,
            base_api=service_config["base_api"],
            base_web=service_config["base_web"],
        )

        print(f"\nAccess Token: {access_token}")
        print(f"Access Token Secret: {access_secret}")

        save_to_credentials_py(
            client_id,
            client_secret,
            access_token,
            access_secret,
            credentials_path,
            credentials_file,
        )

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
