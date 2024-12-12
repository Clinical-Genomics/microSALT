import requests
from microSALT.utils.pubmlst.authentication import generate_oauth_header
from microSALT.utils.pubmlst.helpers import fetch_paginated_data

BASE_API = "https://rest.pubmlst.org"

def validate_session_token(session_token, session_secret):
    """Ensure session token and secret are valid."""
    if not session_token or not session_secret:
        raise ValueError("Session token or secret is missing. Please authenticate first.")

def query_databases(session_token, session_secret):
    """Query available PubMLST databases."""
    validate_session_token(session_token, session_secret)
    url = f"{BASE_API}/db"
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        res = response.json()
        # Ensure the response is a list of database entries
        if not isinstance(res, list):
            raise ValueError(f"Unexpected response format from /db endpoint: {res}")
        return res
    else:
        raise ValueError(f"Failed to query databases: {response.status_code} - {response.text}")


def fetch_schemes(database, session_token, session_secret):
    """Fetch available schemes for a database."""
    validate_session_token(session_token, session_secret)
    url = f"{BASE_API}/db/{database}/schemes"
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError(f"Failed to fetch schemes: {response.status_code} - {response.text}")

def download_profiles(database, scheme_id, session_token, session_secret):
    """Download MLST profiles."""
    validate_session_token(session_token, session_secret)
    if not scheme_id:
        raise ValueError("Scheme ID is required to download profiles.")
    url = f"{BASE_API}/db/{database}/schemes/{scheme_id}/profiles"
    return fetch_paginated_data(url, session_token, session_secret)



def download_locus(database, locus, session_token, session_secret):
    """Download locus sequence files."""
    validate_session_token(session_token, session_secret)
    url = f"{BASE_API}/db/{database}/loci/{locus}/alleles_fasta"
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.content  # Return raw FASTA content
    else:
        raise ValueError(f"Failed to download locus: {response.status_code} - {response.text}")

def check_database_metadata(database, session_token, session_secret):
    """Check database metadata (last update)."""
    validate_session_token(session_token, session_secret)
    url = f"{BASE_API}/db/{database}"
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError(
            f"Failed to check database metadata: {response.status_code} - {response.text}"
        )
