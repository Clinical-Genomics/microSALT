class MicroSALTError(Exception):
    """Base class for exceptions in MicroSALT."""

    def __init__(self, message: str = ""):
        super().__init__(message)


class RefUpdateLockError(MicroSALTError):
    """Raised when an operation is attempted while a reference update is in progress."""
    pass
