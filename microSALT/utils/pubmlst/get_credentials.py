#!/usr/bin/env python3

import json
import os
import sys

from rauth import OAuth1Service

BASE_WEB = {
    "PubMLST": "https://pubmlst.org/bigsdb",
}
BASE_API = {
    "PubMLST": "https://rest.pubmlst.org",
}

SITE = "PubMLST"
DB = "pubmlst_test_seqdef"

# Import client_id and client_secret from credentials.py
try:
    from microSALT.utils.pubmlst_old.credentials import CLIENT_ID, CLIENT_SECRET
except ImportError:
    print("Error: 'credentials.py' file not found or missing CLIENT_ID and CLIENT_SECRET.")
    sys.exit(1)


def main():
    site = SITE
    db = DB

    access_token, access_secret = get_new_access_token(site, db, CLIENT_ID, CLIENT_SECRET)
    print(f"\nAccess Token: {access_token}")
    print(f"Access Token Secret: {access_secret}")

    save_to_credentials_py(CLIENT_ID, CLIENT_SECRET, access_token, access_secret)


def get_new_access_token(site, db, client_id, client_secret):
    """Obtain a new access token and secret."""
    service = OAuth1Service(
        name="BIGSdb_downloader",
        consumer_key=client_id,
        consumer_secret=client_secret,
        request_token_url=f"{BASE_API[site]}/db/{db}/oauth/get_request_token",
        access_token_url=f"{BASE_API[site]}/db/{db}/oauth/get_access_token",
        base_url=BASE_API[site],
    )

    request_token, request_secret = get_request_token(service)
    print(
        "Please log in using your user account at "
        f"{BASE_WEB[site]}?db={db}&page=authorizeClient&oauth_token={request_token} "
        "using a web browser to obtain a verification code."
    )
    verifier = input("Please enter verification code: ")

    # Exchange request token for access token
    raw_access = service.get_raw_access_token(
        request_token, request_secret, params={"oauth_verifier": verifier}
    )
    if raw_access.status_code != 200:
        print(f"Error obtaining access token: {raw_access.text}")
        sys.exit(1)

    access_data = raw_access.json()
    return access_data["oauth_token"], access_data["oauth_token_secret"]


def get_request_token(service):
    """Handle JSON response from the request token endpoint."""
    response = service.get_raw_request_token(params={"oauth_callback": "oob"})
    if response.status_code != 200:
        print(f"Error obtaining request token: {response.text}")
        sys.exit(1)
    try:
        data = json.loads(response.text)
        return data["oauth_token"], data["oauth_token_secret"]
    except json.JSONDecodeError:
        print(f"Failed to parse JSON response: {response.text}")
        sys.exit(1)


def save_to_credentials_py(client_id, client_secret, access_token, access_secret):
    """Save tokens in the credentials.py file."""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    credentials_path = os.path.join(script_dir, "credentials.py")
    with open(credentials_path, "w") as f:
        f.write(f'CLIENT_ID = "{client_id}"\n')
        f.write(f'CLIENT_SECRET = "{client_secret}"\n')
        f.write(f'ACCESS_TOKEN = "{access_token}"\n')
        f.write(f'ACCESS_SECRET = "{access_secret}"\n')
    print(f"Tokens saved to {credentials_path}")


if __name__ == "__main__":
    main()
