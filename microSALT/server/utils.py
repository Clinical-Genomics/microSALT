import base64
from importlib.resources import files
from pathlib import Path


import base64
from importlib.resources import files
from importlib.resources.abc import Traversable

SWEDAC_LOGO_PATH: Traversable = files("microSALT").joinpath("artwork", "swedac.jpg")
MICROSALT_LOGO_PATH: Traversable = files("microSALT").joinpath("artwork", "microsalt.jpg")


def read_jpg(file_path: Traversable) -> str:
    """Return a base64-encoded data URI for a JPEG image."""
    with file_path.open("rb") as jpg_file:
        encoded_string: str = base64.b64encode(jpg_file.read()).decode("utf-8")
    return f"data:image/jpeg;base64,{encoded_string}"
