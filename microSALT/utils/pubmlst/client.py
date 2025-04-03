import requests
from urllib.parse import urlencode
from microSALT.utils.pubmlst.helpers import (
    BASE_API,
    generate_oauth_header,
    load_auth_credentials,
    parse_pubmlst_url,
    should_skip_endpoint,
    get_db_type_capabilities
)
from microSALT.utils.pubmlst.constants import RequestType, HTTPMethod, ResponseHandler
from microSALT.utils.pubmlst.exceptions import PUBMLSTError, SessionTokenRequestError
from microSALT.utils.pubmlst.authentication import load_session_credentials
from microSALT.utils.pubmlst.authentication import get_new_session_token
from microSALT import logger

class PubMLSTClient:
    """Client for interacting with the PubMLST authenticated API."""

    def __init__(self):
        """Initialize the PubMLST client."""
        try:
            self.consumer_key, self.consumer_secret, self.access_token, self.access_secret = load_auth_credentials()
            self.database = "pubmlst_test_seqdef"
            self.session_token, self.session_secret = load_session_credentials(self.database)
        except PUBMLSTError as e:
            logger.error(f"Failed to initialize PubMLST client: {e}")
            raise


    @staticmethod
    def parse_pubmlst_url(url: str):
        """
        Wrapper for the parse_pubmlst_url function.
        """
        return parse_pubmlst_url(url)


    def _make_request(
    self,
    request_type: RequestType,
    method: HTTPMethod,
    url: str,
    db: str = None,
    response_handler: ResponseHandler = ResponseHandler.JSON,
    retry_on_401: bool = True,
    skip_submissions: bool = True
):
        """Handle API requests, support retry on 401, and robust error handling."""
        try:
            if skip_submissions and "/submissions" in url:
                logger.debug(f"[SKIP] Skipping submission-related URL: {url}")
                return None
            
            if db and should_skip_endpoint(url, get_db_type_capabilities(db)):
                logger.debug(f"[SKIP] Skipping incompatible URL for {db}: {url}")
                return None

            if db:
                session_token, session_secret = load_session_credentials(db)
            else:
                session_token, session_secret = self.session_token, self.session_secret

            if request_type == RequestType.AUTH:
                headers = {
                    "Authorization": generate_oauth_header(
                        url,
                        self.consumer_key,
                        self.consumer_secret,
                        self.access_token,
                        self.access_secret,
                    )
                }
            elif request_type == RequestType.DB:
                headers = {
                    "Authorization": generate_oauth_header(
                        url,
                        self.consumer_key,
                        self.consumer_secret,
                        session_token,
                        session_secret,
                    )
                }
            else:
                raise ValueError(f"Unsupported request type: {request_type}")


            response = requests.request(method.value, url, headers=headers)

            if response.status_code == 401 and retry_on_401 and db:
                logger.debug(f"[DEBUG] Got 401 Unauthorized. Refreshing session token and retrying for {url}")
                get_new_session_token(db)
                return self._make_request(
                    request_type=request_type,
                    method=method,
                    url=url,
                    db=db,
                    response_handler=response_handler,
                    retry_on_401=False 
                )

            if response.status_code == 404:
                logger.debug(f"[DEBUG] 404 Not Found for {url}")
                return None

            response.raise_for_status()

            if response_handler == ResponseHandler.CONTENT:
                return response.content
            elif response_handler == ResponseHandler.TEXT:
                return response.text
            elif response_handler == ResponseHandler.JSON:
                return response.json()
            else:
                raise ValueError(f"Unsupported response handler: {response_handler}")

        except requests.exceptions.HTTPError as e:
            raise SessionTokenRequestError(db or self.database, e.response.status_code, e.response.text) from e
        except requests.exceptions.RequestException as e:
            logger.error(f"Request failed: {e}")
            raise PUBMLSTError(f"Request failed: {e}") from e
        except Exception as e:
            logger.error(f"Unexpected error during request: {e}")
            raise PUBMLSTError(f"An unexpected error occurred: {e}") from e


    def query_databases(self):
        """Query available PubMLST databases."""
        url = f"{BASE_API}/db"
        return self._make_request(RequestType.DB, HTTPMethod.GET, url, response_handler=ResponseHandler.JSON)


    def download_locus(self, db: str, locus: str, **kwargs):
        """Download locus sequence files."""
        base_url = f"{BASE_API}/db/{db}/loci/{locus}/alleles_fasta"
        query_string = urlencode(kwargs)
        url = f"{base_url}?{query_string}" if query_string else base_url
        return self._make_request(RequestType.DB, HTTPMethod.GET, url, db=db, response_handler=ResponseHandler.TEXT)


    def download_profiles_csv(self, db: str, scheme_id: int):
        """Download MLST profiles in CSV format."""
        if not scheme_id:
            raise ValueError("Scheme ID is required to download profiles CSV.")
        url = f"{BASE_API}/db/{db}/schemes/{scheme_id}/profiles_csv"
        return self._make_request(RequestType.DB, HTTPMethod.GET, url, db=db, response_handler=ResponseHandler.TEXT)


    def retrieve_scheme_info(self, db: str, scheme_id: int):
        """Retrieve information about a specific MLST scheme."""
        url = f"{BASE_API}/db/{db}/schemes/{scheme_id}"
        return self._make_request(RequestType.DB, HTTPMethod.GET, url, db=db, response_handler=ResponseHandler.JSON)


    def list_schemes(self, db: str):
        """List available MLST schemes for a specific database."""
        url = f"{BASE_API}/db/{db}/schemes"
        return self._make_request(RequestType.DB, HTTPMethod.GET, url, db=db, response_handler=ResponseHandler.JSON)
