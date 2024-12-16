from enum import Enum

class RequestType(Enum):
    AUTH = "auth"
    DB = "db"

class CredentialsFile(Enum):
    MAIN = "main"
    SESSION = "session"

class Encoding(Enum):
    UTF8 = "utf-8"

class HTTPMethod(Enum):
    GET = "GET"
    POST = "POST"
    PUT = "PUT"
    DELETE = "DELETE"
    PATCH = "PATCH"
    HEAD = "HEAD"
    OPTIONS = "OPTIONS"

class ResponseHandler(Enum):
    CONTENT = "content"
    TEXT = "text"
    JSON = "json"