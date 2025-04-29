import sys
import os

from argparse import ArgumentParser
from rauth import OAuth1Service
from microSALT import app
from microSALT.utils.pubmlst.helpers import get_path, get_service_config, folders_config
from microSALT.utils.pubmlst.constants import CREDENTIALS_KEY

db = "pubmlst_test_seqdef"


def validate_credentials(client_id, client_secret):
    """Ensure client_id and client_secret are not empty."""
    if not client_id or not client_id.strip():
        raise ValueError("Invalid CLIENT_ID: It must not be empty.")
    if not client_secret or not client_secret.strip():
        raise ValueError("Invalid CLIENT_SECRET: It must not be empty.")


def get_request_token(service):
    """Handle JSON response from the request token endpoint."""
    response = service.get_raw_request_token(params={"oauth_callback": "oob"})
    if not response.ok:
        print(f"Error obtaining request token: {response.text}")
        sys.exit(1)
    data = response.json()
    return data["oauth_token"], data["oauth_token_secret"]


def get_new_access_token(client_id, client_secret, db: str, base_api: str, base_web: str):
    """Obtain a new access token and secret."""
    service = OAuth1Service(
        name="BIGSdb_downloader",
        consumer_key=client_id,
        consumer_secret=client_secret,
        request_token_url=f"{base_api}/db/{db}/oauth/get_request_token",
        access_token_url=f"{base_api}/db/{db}/oauth/get_access_token",
        base_url=base_api,
    )
    request_token, request_secret = get_request_token(service)
    print(
        "Please log in using your user account at "
        f"{base_web}?db={db}&page=authorizeClient&oauth_token={request_token} "
        "using a web browser to obtain a verification code."
    )
    verifier = input("Please enter verification code: ")

    raw_access = service.get_raw_access_token(
        request_token, request_secret, params={"oauth_verifier": verifier}
    )
    if not raw_access.ok:
        print(f"Error obtaining access token: {raw_access.text}")
        sys.exit(1)

    access_data = raw_access.json()
    return access_data["oauth_token"], access_data["oauth_token_secret"]


def save_to_credentials_py(
    client_id, client_secret, access_token, access_secret, credentials_path, credentials_file
):
    """Save tokens in the credentials.py file."""
    credentials_path.mkdir(parents=True, exist_ok=True)

    with open(credentials_file, "w") as f:
        f.write(f'CLIENT_ID = "{client_id}"\n')
        f.write(f'CLIENT_SECRET = "{client_secret}"\n')
        f.write(f'ACCESS_TOKEN = "{access_token}"\n')
        f.write(f'ACCESS_SECRET = "{access_secret}"\n')
    print(f"Tokens saved to {credentials_file}")


def main(service, species=None):
    try:
        service_config = get_service_config(service)
        bigsd_config = service_config["config"]
        client_id = bigsd_config["client_id"]
        client_secret = bigsd_config["client_secret"]
        validate_credentials(client_id, client_secret)

        # Determine the database
        if service == "pubmlst":
            db = service_config["database"]
        elif service == "pasteur":
            if not species:
                raise ValueError("For the 'pasteur' service, you must provide a species.")
            db = f"pasteur_{species}_seqdef"
        else:
            raise ValueError(f"Unknown service: {service}")

        credentials_path = get_path(folders_config, CREDENTIALS_KEY)
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


if __name__ == "__main__":
    parser = ArgumentParser(description="Get PUBMLST or Pasteur credentials.")
    parser.add_argument(
        "-s",
        "--service",
        type=str,
        default="pubmlst",
        help="Service name (default: pubmlst)",
    )
    parser.add_argument(
        "-sp",
        "--species",
        type=str,
        help="Species name (required for the 'pasteur' service)",
    )
    args = parser.parse_args()
    main(args.service, args.species)
