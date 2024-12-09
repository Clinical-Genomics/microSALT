import requests

from microSALT.utils.pubmlst.authentication import generate_oauth_header


def fetch_paginated_data(url, session_token, session_secret):
    """Fetch paginated data using the session token and secret."""
    results = []
    while url:
        headers = {"Authorization": generate_oauth_header(url, session_token, session_secret)}
        response = requests.get(url, headers=headers)

        print(f"Fetching URL: {url}")
        print(f"Response Status Code: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            results.extend(data.get("profiles", []))
            url = data.get("paging", {}).get("next", None)  # Get the next page URL if available
        else:
            raise ValueError(f"Failed to fetch data: {response.status_code} - {response.text}")
    return results
