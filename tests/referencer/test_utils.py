import pytest
from microSALT.utils.referencer.utils import (
    get_species,
    is_escherichia,
    is_klebsiella,
    get_reference_if_enterobacteriaceae,
)
from microSALT.utils.referencer.constants import ESCHERICHIA_REFERENCE, KLEBSIELLA_REFERENCE, Escherichia, Klebsiella

@pytest.mark.parametrize(
    "organism_name, expected_species",
    [
        ("Escherichia coli", "coli"),
        ("Klebsiella pneumoniae", "pneumoniae"),
        ("Escherichia coli O157:H7", "o157:h7"),
    ],
)
def test_get_species(organism_name, expected_species):
    assert get_species(organism_name) == expected_species

@pytest.mark.parametrize(
    "species_name, expected_result",
    [
        (Escherichia.ALBA, True),
        (Escherichia.COLI, True),
        (Klebsiella.ELECTRICA, False),
    ],
)
def test_is_escherichia(species_name, expected_result):
    assert is_escherichia(species_name) == expected_result

@pytest.mark.parametrize(
    "species_name, expected_result",
    [
        (Klebsiella.PNEUMONIAE, True),
        (Klebsiella.OXYTOCA, True),
        (Escherichia.ALBA, False),
    ],
)
def test_is_klebsiella(species_name, expected_result):
    assert is_klebsiella(species_name) == expected_result

@pytest.mark.parametrize(
    "organism_name, expected_reference",
    [
        (Escherichia.COLI, ESCHERICHIA_REFERENCE),
        (Klebsiella.PNEUMONIAE, KLEBSIELLA_REFERENCE),
        ("other species", "other species"),
    ],
)
def test_get_reference_if_enterobacteriaceae(organism_name, expected_reference):
    assert get_reference_if_enterobacteriaceae(organism_name) == expected_reference