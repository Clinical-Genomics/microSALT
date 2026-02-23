class MicroSALTError(Exception):
    """Base class for exceptions in MicroSALT."""

    def __init__(self, message: str = ""):
        super().__init__(message)
