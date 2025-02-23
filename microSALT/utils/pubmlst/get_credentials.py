import sys
import os
from rauth import OAuth1Service
from microSALT import app
from microSALT.utils.pubmlst.helpers import get_path, BASE_API, BASE_WEB, folders_config, credentials_path_key, pubmlst_auth_credentials_file_name

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


def get_new_access_token(client_id, client_secret):
    """Obtain a new access token and secret."""
    service = OAuth1Service(
        name="BIGSdb_downloader",
        consumer_key=client_id,
        consumer_secret=client_secret,
        request_token_url=f"{BASE_API}/db/{db}/oauth/get_request_token",
        access_token_url=f"{BASE_API}/db/{db}/oauth/get_access_token",
        base_url=BASE_API,
    )
    request_token, request_secret = get_request_token(service)
    print(
        "Please log in using your user account at "
        f"{BASE_WEB}?db={db}&page=authorizeClient&oauth_token={request_token} "
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


def save_to_credentials_py(client_id, client_secret, access_token, access_secret, credentials_path, credentials_file):
    """Save tokens in the credentials.py file."""
    credentials_path.mkdir(parents=True, exist_ok=True)

    with open(credentials_file, "w") as f:
        f.write(f'CLIENT_ID = "{client_id}"\n')
        f.write(f'CLIENT_SECRET = "{client_secret}"\n')
        f.write(f'ACCESS_TOKEN = "{access_token}"\n')
        f.write(f'ACCESS_SECRET = "{access_secret}"\n')
    print(f"Tokens saved to {credentials_file}")


def main():
    try:
        pubmlst_config = app.config["pubmlst"]
        client_id = pubmlst_config["client_id"]
        client_secret = pubmlst_config["client_secret"]
        validate_credentials(client_id, client_secret)
        credentials_path = get_path(folders_config, credentials_path_key)
        credentials_file = os.path.join(credentials_path, pubmlst_auth_credentials_file_name)
        access_token, access_secret = get_new_access_token(client_id, client_secret)
        print(f"\nAccess Token: {access_token}")
        print(f"Access Token Secret: {access_secret}")
        save_to_credentials_py(client_id, client_secret, access_token, access_secret, credentials_path, credentials_file)

    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
