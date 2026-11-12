class MicroSALTError(Exception):
    """Base class for exceptions in MicroSALT."""

    def __init__(self, message: str = ""):
        super().__init__(message)


class RefUpdateLockError(MicroSALTError):
    """Raised when an operation is attempted while a reference update is in progress."""

    pass


class JobCreationError(MicroSALTError):
    """Raised when there is an error creating a job."""

    pass


class ProfileCreationError(MicroSALTError):
    """Raised when there is an error creating a profile."""

    pass


class ReferenceSetupError(MicroSALTError):
    """A required reference could not be fetched, validated or indexed."""


class AnalysisSetupError(MicroSALTError):
    """Collected setup failures that must prevent SLURM submission."""

    def __init__(self, errors: list[str]):
        self.errors = errors
        super().__init__(
            "Analysis setup failed. No SLURM jobs were submitted.\n\n"
            + "\n".join(f"- {error}" for error in errors)
        )
