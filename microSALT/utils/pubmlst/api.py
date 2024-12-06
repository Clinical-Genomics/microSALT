import requests

from microSALT.utils.pubmlst.authentication import generate_oauth_header
from microSALT.utils.pubmlst.helpers import fetch_paginated_data

BASE_API = "https://rest.pubmlst.org"


def query_databases(session_token, session_secret):
    """Query available PubMLST databases."""
    url = f"{BASE_API}/db"
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError(f"Failed to query databases: {response.status_code} - {response.text}")


def fetch_schemes(database, session_token, session_secret):
    """Fetch available schemes for a database."""
    url = f"{BASE_API}/db/{database}/schemes"
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError(f"Failed to fetch schemes: {response.status_code} - {response.text}")


def download_profiles(database, scheme_id, session_token, session_secret):
    """Download MLST profiles."""
    url = f"{BASE_API}/db/{database}/schemes/{scheme_id}/profiles"
    return fetch_paginated_data(url, session_token, session_secret)


def download_locus(database, locus, session_token, session_secret):
    """Download locus sequence files."""
    url = f"{BASE_API}/db/{database}/loci/{locus}/alleles_fasta"
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.content  # Return raw FASTA content
    else:
        raise ValueError(f"Failed to download locus: {response.status_code} - {response.text}")


def check_database_metadata(database, session_token, session_secret):
    """Check database metadata (last update)."""
    url = f"{BASE_API}/db/{database}"
    headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
    response = requests.get(url, headers=headers)
    if response.status_code == 200:
        return response.json()
    else:
        raise ValueError(
            f"Failed to check database metadata: {response.status_code} - {response.text}"
        )
