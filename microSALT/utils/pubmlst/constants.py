from enum import Enum
from werkzeug.routing import Map, Rule

def add_prefix_to_rules(url_map: Map, prefix: str):
    """
    Add a prefix to all rules in the given URL map.

    :param url_map: The original URL map.
    :param prefix: The prefix to add (e.g., '/api').
    :return: A new URL map with the prefix applied.
    """
    new_rules = []
    new_rules.extend(
        Rule(
            f"{prefix}{rule.rule}",
            endpoint=rule.endpoint,
            methods=rule.methods,
        )
        for rule in url_map.iter_rules()
    )
    return Map(new_rules)

CREDENTIALS_KEY: str = "credentials"


class RequestType(Enum):
    AUTH = "auth"
    DB = "db"


class CredentialsFile(Enum):
    MAIN = "main"
    SESSION = "session"


class Encoding(Enum):
    UTF8 = "utf-8"


class HTTPMethod(Enum):
    GET = "GET"
    POST = "POST"
    PUT = "PUT"
    DELETE = "DELETE"
    PATCH = "PATCH"
    HEAD = "HEAD"
    OPTIONS = "OPTIONS"


class ResponseHandler(Enum):
    CONTENT = "content"
    TEXT = "text"
    JSON = "json"


pubmlst_urls = Map(
    [
        Rule("/", endpoint="root"),
        Rule("/db", endpoint="db_root"),
        Rule("/db/<db>", endpoint="database_root"),
        Rule("/db/<db>/classification_schemes", endpoint="classification_schemes"),
        Rule(
            "/db/<db>/classification_schemes/<int:classification_scheme_id>",
            endpoint="classification_scheme",
        ),
        Rule(
            "/db/<db>/classification_schemes/<int:classification_scheme_id>/groups",
            endpoint="classification_scheme_groups",
        ),
        Rule(
            "/db/<db>/classification_schemes/<int:classification_scheme_id>/groups/<int:group_id>",
            endpoint="classification_scheme_group",
        ),
        Rule("/db/<db>/loci", endpoint="loci"),
        Rule("/db/<db>/loci/<locus>", endpoint="locus"),
        Rule("/db/<db>/loci/<locus>/alleles", endpoint="locus_alleles"),
        Rule("/db/<db>/loci/<locus>/alleles_fasta", endpoint="locus_alleles_fasta"),
        Rule("/db/<db>/loci/<locus>/alleles/<int:allele_id>", endpoint="locus_allele"),
        Rule("/db/<db>/loci/<locus>/sequence", endpoint="locus_sequence_post"),
        Rule("/db/<db>/sequence", endpoint="sequence_post"),
        Rule("/db/<db>/sequences", endpoint="sequences"),
        Rule("/db/<db>/schemes", endpoint="schemes"),
        Rule("/db/<db>/schemes/<int:scheme_id>", endpoint="scheme"),
        Rule("/db/<db>/schemes/<int:scheme_id>/loci", endpoint="scheme_loci"),
        Rule("/db/<db>/schemes/<int:scheme_id>/fields/<field>", endpoint="scheme_field"),
        Rule("/db/<db>/schemes/<int:scheme_id>/profiles", endpoint="scheme_profiles"),
        Rule("/db/<db>/schemes/<int:scheme_id>/profiles_csv", endpoint="scheme_profiles_csv"),
        Rule(
            "/db/<db>/schemes/<int:scheme_id>/profiles/<int:profile_id>", endpoint="scheme_profile"
        ),
        Rule("/db/<db>/schemes/<int:scheme_id>/sequence", endpoint="scheme_sequence_post"),
        Rule("/db/<db>/schemes/<int:scheme_id>/designations", endpoint="scheme_designations_post"),
        Rule("/db/<db>/isolates", endpoint="isolates"),
        Rule("/db/<db>/genomes", endpoint="genomes"),
        Rule("/db/<db>/isolates/search", endpoint="isolates_search_post"),
        Rule("/db/<db>/isolates/<int:isolate_id>", endpoint="isolate"),
        Rule(
            "/db/<db>/isolates/<int:isolate_id>/allele_designations",
            endpoint="isolate_allele_designations",
        ),
        Rule(
            "/db/<db>/isolates/<int:isolate_id>/allele_designations/<locus>",
            endpoint="isolate_allele_designation_locus",
        ),
        Rule("/db/<db>/isolates/<int:isolate_id>/allele_ids", endpoint="isolate_allele_ids"),
        Rule(
            "/db/<db>/isolates/<int:isolate_id>/schemes/<int:scheme_id>/allele_designations",
            endpoint="isolate_scheme_allele_designations",
        ),
        Rule(
            "/db/<db>/isolates/<int:isolate_id>/schemes/<int:scheme_id>/allele_ids",
            endpoint="isolate_scheme_allele_ids",
        ),
        Rule("/db/<db>/isolates/<int:isolate_id>/contigs", endpoint="isolate_contigs"),
        Rule("/db/<db>/isolates/<int:isolate_id>/contigs_fasta", endpoint="isolate_contigs_fasta"),
        Rule("/db/<db>/isolates/<int:isolate_id>/history", endpoint="isolate_history"),
        Rule("/db/<db>/contigs/<int:contig_id>", endpoint="contig"),
        Rule("/db/<db>/fields", endpoint="fields"),
        Rule("/db/<db>/fields/<field>", endpoint="field"),
        Rule("/db/<db>/users/<int:user_id>", endpoint="user"),
        Rule("/db/<db>/curators", endpoint="curators"),
        Rule("/db/<db>/projects", endpoint="projects"),
        Rule("/db/<db>/projects/<int:project_id>", endpoint="project"),
        Rule("/db/<db>/projects/<int:project_id>/isolates", endpoint="project_isolates"),
        Rule("/db/<db>/submissions", endpoint="submissions"),
        Rule("/db/<db>/submissions/<int:submission_id>", endpoint="submission"),
        Rule("/db/<db>/submissions/<int:submission_id>/messages", endpoint="submission_messages"),
        Rule("/db/<db>/submissions/<int:submission_id>/files", endpoint="submission_files"),
        Rule(
            "/db/<db>/submissions/<int:submission_id>/files/<filename>", endpoint="submission_file"
        ),
    ]
)

pasteur_urls = add_prefix_to_rules(url_map, "/api")

URL_MAPS = {
    "pubmlst": pubmlst_urls,
    "pasteur": pasteur_urls,
}
