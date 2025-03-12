from microSALT.utils.referencer.constants import Escherichia, Klebsiella


def is_escherichia(organism_name: str) -> bool:
    """Checks if the organism name is a proper Escherichia."""
    organism_name = organism_name.lower()
    return any(organism_name == member.value for member in Escherichia)


def is_klebsiella(organism_name: str) -> bool:
    """Checks if the organism name is a proper Klebsiella."""
    organism_name = organism_name.lower()
    return any(organism_name == member.value for member in Klebsiella)


def get_reference_if_enterobacteriaceae(organism_name: str) -> str:
    """Returns the reference for the organism name if it belongs to Enterobacteriaceae."""
    GENUS_REFERENCE_MAP = {
        "escherichia": (is_escherichia, Escherichia.ESCHERICHIA_REFERENCE),
        "klebsiella": (is_klebsiella, Klebsiella.KLEBSIELLA_REFERENCE),
    }
    organism_name = organism_name.lower()
    return next(
        (
            reference
            for check_func, reference in GENUS_REFERENCE_MAP.values()
            if check_func(organism_name)
        ),
        organism_name,
    )
