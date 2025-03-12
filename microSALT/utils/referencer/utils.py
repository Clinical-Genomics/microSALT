from microSALT.utils.referencer.constants import (
    Escherichia,
    Klebsiella,
    ESCHERICHIA_REFERENCE,
    KLEBSIELLA_REFERENCE,
)


def get_species(organism_name: str) -> str:
    """Returns the species of the organism name."""
    return organism_name.lower().split()[-1]


def is_escherichia(species_name: str) -> bool:
    """Checks if the organism name is a proper Escherichia."""
    return any(species_name == member.value for member in Escherichia)


def is_klebsiella(species_name: str) -> bool:
    """Checks if the organism name is a proper Klebsiella."""
    return any(species_name == member.value for member in Klebsiella)


def get_reference_if_enterobacteriaceae(organism_name: str) -> str:
    """Returns the reference for the organism name if it belongs to Enterobacteriaceae."""
    species: str = get_species(organism_name)
    GENUS_REFERENCE_MAP = {
        "escherichia": (is_escherichia, ESCHERICHIA_REFERENCE),
        "klebsiella": (is_klebsiella, KLEBSIELLA_REFERENCE),
    }
    return next(
        (
            reference
            for check_func, reference in GENUS_REFERENCE_MAP.values()
            if check_func(species_name=species)
        ),
        organism_name,
    )
