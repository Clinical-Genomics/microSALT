"""Suite registry. Release-specific suites live in sub-packages named after the release."""

from importlib import import_module
from pkgutil import iter_modules

from core.suite import Suite

from .stability import StabilitySuite

REGISTRY: dict[str, type[Suite]] = {StabilitySuite.name: StabilitySuite}


def _register_package(package: str) -> None:
    module = import_module(package)
    for info in iter_modules(module.__path__):
        submodule = import_module(f"{package}.{info.name}")
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
