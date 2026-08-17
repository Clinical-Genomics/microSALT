import base64
from importlib.resources import files
from pathlib import Path


def get_project_root_dir() -> Path:
    return Path(files("microSALT"))


SWEDAC_LOGO_PATH = Path(
    get_project_root_dir(),
    "artwork",
    "swedac.jpg",
)
MICROSALT_LOGO_PATH = Path(
    get_project_root_dir(),
    "artwork",
    "microsalt.jpg",
)


def read_jpg(file_path: Path) -> str:
    """Return base64 encoding of a PNG image."""
    with open(file_path, "rb") as jpg_file:
        encoded_string: str = base64.b64encode(jpg_file.read()).decode("utf-8")
    return f"data:image/jpeg;base64,{encoded_string}"
