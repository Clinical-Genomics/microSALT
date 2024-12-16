class PUBMLSTError(Exception):
    """Base exception for PUBMLST utilities."""
    def __init__(self, message=None):
        super(PUBMLSTError, self).__init__(f"PUBMLST: {message}")


class CredentialsFileNotFound(PUBMLSTError):
    """Raised when the PUBMLST credentials file is not found."""
    def __init__(self, credentials_file):
        message = (
            f"Credentials file not found: {credentials_file}. "
            "Please generate it using the get_credentials script."
        )
        super(CredentialsFileNotFound, self).__init__(message)


class InvalidCredentials(PUBMLSTError):
    """Raised when the credentials file contains invalid or missing fields."""
    def __init__(self, missing_fields):
        message = (
            "Invalid credentials: All fields (CLIENT_ID, CLIENT_SECRET, ACCESS_TOKEN, ACCESS_SECRET) "
            f"must be non-empty. Missing or empty fields: {', '.join(missing_fields)}. "
            "Please regenerate the credentials file using the get_credentials script."
        )
        super(InvalidCredentials, self).__init__(message)


class PathResolutionError(PUBMLSTError):
    """Raised when the file path cannot be resolved from the configuration."""
    def __init__(self, config_key):
        message = (
            f"Failed to resolve the path for configuration key: '{config_key}'. "
            "Ensure it is correctly set in the configuration."
        )
        super(PathResolutionError, self).__init__(message)


class SaveSessionError(PUBMLSTError):
    """Raised when saving the session token fails."""
    def __init__(self, db, reason):
        message = f"Failed to save session token for database '{db}': {reason}"
        super(SaveSessionError, self).__init__(message)


class SessionTokenRequestError(PUBMLSTError):
    """Raised when requesting a session token fails."""
    def __init__(self, db, status_code, response_text):
        message = f"Failed to fetch session token for database '{db}': {status_code} - {response_text}"
        super(SessionTokenRequestError, self).__init__(message)


class SessionTokenResponseError(PUBMLSTError):
    """Raised when the session token response is invalid."""
    def __init__(self, db, reason):
        message = f"Invalid session token response for database '{db}': {reason}"
        super(SessionTokenResponseError, self).__init__(message)
