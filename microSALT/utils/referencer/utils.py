from microSALT.utils.referencer.constants import (
    Escherichia,
    Klebsiella,
    ESCHERICHIA_REFERENCE,
    KLEBSIELLA_REFERENCE,
)


def get_species(organism_name: str) -> str:
    """Returns the species of the organism name."""
    return organism_name.lower().split()[-1]


def get_genus(organism_name: str) -> str:
    """Returns the genus of the organism name."""
    return organism_name.lower().split()[0]


def is_enterobacteriaceae(organism_name: str) -> bool:
    """Checks if the organism name belongs to Enterobacteriaceae."""
    return get_genus(organism_name) in ("escherichia", "klebsiella")


def is_escherichia(species_name) -> bool:
    """Checks if the organism name is a proper Escherichia."""
    return any(species_name == member.value for member in Escherichia)


def is_klebsiella(species_name) -> bool:
    """Checks if the organism name is a proper Klebsiella."""
    return any(species_name == member.value for member in Klebsiella)


def get_reference_if_enterobacteriaceae(organism_name: str) -> str:
    """Returns the reference for the organism name if it belongs to Enterobacteriaceae."""
    GENUS_REFERENCE_MAP = {
        (is_escherichia, ESCHERICHIA_REFERENCE),
        (is_klebsiella, KLEBSIELLA_REFERENCE),
    }

    species: str = get_species(organism_name)
    if is_enterobacteriaceae(organism_name):
        for is_species, reference in GENUS_REFERENCE_MAP:
            if is_species(species):
                return reference
    return organism_name
