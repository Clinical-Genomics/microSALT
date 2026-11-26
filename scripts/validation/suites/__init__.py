"""Suite registry. Release-specific suites live in sub-packages named after the release.

A suite module that fails to import is recorded in BROKEN instead of taking down the registry,
so `run_validation.py suites` can still report what is available and why the rest is not.
"""

import traceback
from importlib import import_module
from pkgutil import iter_modules

from core.suite import Suite

from .stability import StabilitySuite

REGISTRY: dict[str, type[Suite]] = {StabilitySuite.name: StabilitySuite}
BROKEN: dict[str, str] = {}


def _register_package(package: str) -> None:
    module = import_module(package)
    for info in iter_modules(module.__path__):
        qualified = f"{package}.{info.name}"
        try:
            submodule = import_module(qualified)
        except (
            Exception
        ) as error:  # noqa: BLE001 - any import failure must be reported, not raised
            BROKEN[qualified] = "".join(
                traceback.format_exception_only(type(error), error)
            ).strip()
            continue
        for value in vars(submodule).values():
            if (
                isinstance(value, type)
                and issubclass(value, Suite)
                and value is not Suite
                and hasattr(value, "name")
            ):
                REGISTRY[value.name] = value


for _info in iter_modules(__path__):
    if _info.ispkg:
        _register_package(f"{__name__}.{_info.name}")
