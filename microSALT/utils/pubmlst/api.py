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
    url = "{}/db".format(BASE_API)
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        res = response.json()
        if not isinstance(res, list):
            raise ValueError("Unexpected response format from /db endpoint: {}".format(res))
        return res
    else:
        raise ValueError("Failed to query databases: {} - {}".format(response.status_code, response.text))

def fetch_schemes(database, session_token, session_secret):
    """Fetch available schemes for a database."""
    validate_session_token(session_token, session_secret)
    url = "{}/db/{}/schemes".format(BASE_API, database)
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError("Failed to fetch schemes: {} - {}".format(response.status_code, response.text))

def download_profiles_csv(database, scheme_id, session_token, session_secret):
    """Download MLST profiles in CSV format."""
    validate_session_token(session_token, session_secret)
    if not scheme_id:
        raise ValueError("Scheme ID is required to download profiles CSV.")
    url = "{}/db/{}/schemes/{}/profiles_csv".format(BASE_API, database, scheme_id)
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        raise ValueError("Failed to download profiles CSV: {}".format(e))

def download_locus(database, locus, session_token, session_secret):
    """Download locus sequence files."""
    validate_session_token(session_token, session_secret)
    url = "{}/db/{}/loci/{}/alleles_fasta".format(BASE_API, database, locus)
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        return response.content
    except requests.exceptions.RequestException as e:
        raise ValueError("Failed to download locus: {}".format(e))

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