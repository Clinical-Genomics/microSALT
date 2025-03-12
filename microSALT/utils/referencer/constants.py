from enum import Enum


class StrEnum(str, Enum):
    def __str__(self):
        return self.value


class Escherichia(StrEnum):
    ALBA: str = "alba"
    ALBERTII: str = "albertii"
    COLI: str = "coli"
    FAECALIS: str = "faecalis"
    FERGUSONII: str = "fergusonii"
    HOMINIS: str = "hominis"
    MARMOTAE: str = "marmotae"
    RUYSIAE: str = "ruysiae"
    SENEGALENSIS: str = "senegalensis"
    SPP: str = "spp."
    WHITTAMII: str = "whittamii"



class Klebsiella(StrEnum):
    AEROGENES: str = "aerogenes"
    AFRICANA: str = "africana"
    ELECTRICA: str = "electrica"
    GRIMONTII: str = "grimontii"
    HUAXIENSIS: str = "huaxiensis"
    INDICA: str = "indica"
    KIELENSIS: str = "kielensis"
    MICHIGANENSIS: str = "michiganensis"
    MILLETIS: str = "milletis"
    OXYTOCA: str = "oxytoca"
    PASTEURII: str = "pasteurii"
    PNEUMONIAE: str = "pneumoniae"
    OZAENAE: str = "ozaenae"
    GRANULOMATIS: str = "granulomatis"
    SENEGALENSIS: str = "senegalensis"
    SPALLANZANII: str = "spallanzanii"
    SPP: str = "spp."


ESCHERICHIA_REFERENCE = "escherichia_spp."
KLEBSIELLA_REFERENCE = "klebsiella_spp."