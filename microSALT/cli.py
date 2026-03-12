"""This is the main entry point of microSALT.
By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import json
import logging
import os
import pathlib
import re
import subprocess
import sys

import click

from microSALT import __version__, logging_levels, setup_logger
from microSALT.config import MicroSALTConfig, load_config
from microSALT.exc.exceptions import RefUpdateLockError
from microSALT.store.database import (
    create_tables,
    get_scoped_session_registry,
    initialize_database,
)
from microSALT.utils.job_creator import Job_Creator
from microSALT.utils.pubmlst.get_credentials import main as get_bigsdb_credentials_main
from microSALT.utils.referencer import Referencer
from microSALT.utils.reporter import Reporter
from microSALT.utils.scraper import Scraper

default_sampleinfo = {
    "CG_ID_project": "XXX0000",
    "CG_ID_sample": "XXX0000A1",
    "Customer_ID_project": "100100",
    "Customer_ID_sample": "10XY123456",
    "Customer_ID": "cust000",
    "application_tag": "SOMTIN100",
    "date_arrival": "0001-01-01 00:00:00",
    "date_libprep": "0001-01-01 00:00:00",
    "date_sequencing": "0001-01-01 00:00:00",
    "method_libprep": "Not in LIMS",
    "method_sequencing": "Not in LIMS",
    "organism": "Staphylococcus aureus",
    "priority": "standard",
    "reference": "None",
}

logger = logging.getLogger("main_logger")


def done():
    click.echo("INFO - Execution finished!")
    logger.debug("INFO - Execution finished!")


def review_sampleinfo(pfile):
    """Reviews sample info. Returns loaded json object"""

    try:
        with open(pfile) as json_file:
            data = json.load(json_file)
    except Exception:
        click.echo("Unable to read provided sample info file as json. Exiting..")
        sys.exit(-1)

    if isinstance(data, list):
        for entry in data:
            for k, v in default_sampleinfo.items():
                if k not in entry:
                    click.echo(
                        f"WARNING - Parameter {k} needs to be provided in sample json. Formatting example: ({v})"
                    )
    else:
        for k, v in default_sampleinfo.items():
            if k not in data:
                click.echo(
                    f"WARNING - Parameter {k} needs to be provided in sample json. Formatting example: ({v})"
                )
    return data


def teardown_session():
    """Ensure that the session is closed and all resources are released to the connection pool."""
    registry = get_scoped_session_registry()
    if registry:
        registry.remove()


def _ensure_directories(config: MicroSALTConfig) -> None:
    """Create any configured directory paths that do not yet exist."""
    db_uri = config.database.SQLALCHEMY_DATABASE_URI
    db_match = re.search("sqlite:///(.+)", db_uri)
    if db_match:
        db_file = db_match.group(1)
        db_dir = os.path.dirname(db_file)
        if db_dir and not pathlib.Path(db_dir).exists():
            os.makedirs(db_dir)
        proc = subprocess.Popen(f"touch {db_file}".split(), stdout=subprocess.PIPE)
        _, _ = proc.communicate()
        if proc.returncode != 0:
            click.echo("ERROR - Database writing failed! Invalid user access detected!")
            sys.exit(-1)

    folder_paths = [
        config.folders.results,
        config.folders.reports,
        config.folders.seqdata,
        config.folders.profiles,
        config.folders.references,
        config.folders.resistances,
        config.folders.genomes,
        config.folders.credentials,
    ]
    for path in folder_paths:
        p = pathlib.Path(os.path.expandvars(os.path.expanduser(path)))
        if not p.exists():
            os.makedirs(p)


pass_config = click.make_pass_decorator(MicroSALTConfig, ensure=True)


@click.group()
@click.version_option(__version__)
@click.option(
    "--config",
    required=True,
    help="Path to microSALT config JSON file",
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--logging-level",
    default="INFO",
    type=click.Choice(list(logging_levels.keys())),
    help="Set the logging level for the CLI",
)
@click.pass_context
def root(ctx, config, logging_level):
    """microbial Sequence Analysis and Loci-based Typing (microSALT) pipeline"""
    cfg = load_config(config)
    initialize_database(cfg.database.SQLALCHEMY_DATABASE_URI)
    setup_logger(logging_level=logging_level)
    logger.setLevel(logging_levels[logging_level])
    for handler in logger.handlers:
        handler.setLevel(logging_levels[logging_level])
    ctx.obj = cfg
    ctx.call_on_close(teardown_session)


@root.command()
@pass_config
def setup(config: MicroSALTConfig):
    """Create all configured directories and verify database access. Run once after installation."""
    _ensure_directories(config)
    click.echo("INFO - Directory setup complete.")
    create_tables()
    click.echo("INFO - Database tables created (or already exist).")
    done()


@root.command()
@click.argument("sampleinfo_file")
@click.option("--input", help="Full path to input folder", default="")
@click.option(
    "--dry",
    help="Builds instance without posting to SLURM",
    default=False,
    is_flag=True,
)
@click.option("--email", default="", help="Forced e-mail recipient")
@click.option("--skip_update", default=False, help="Skips downloading of references", is_flag=True)
@click.option(
    "--force_update",
    default=False,
    help="Forces downloading of pubMLST references",
    is_flag=True,
)
@click.option("--untrimmed", help="Use untrimmed input data", default=False, is_flag=True)
@pass_config
def analyse(
    config: MicroSALTConfig,
    sampleinfo_file,
    input,
    dry,
    email,
    skip_update,
    force_update,
    untrimmed,
):
    """Sequence analysis, typing and resistance identification"""
    pool = []
    if email:
        config.regex.mail_recipient = email
    config.dry = dry
    if not os.path.isdir(input):
        click.echo(f"ERROR - Sequence data folder {input} does not exist.")
        click.Abort()
    for subfolder in os.listdir(input):
        if os.path.isdir(f"{input}/{subfolder}"):
            pool.append(subfolder)

    run_settings = {
        "input": input,
        "dry": dry,
        "email": config.regex.mail_recipient,
        "skip_update": skip_update,
        "trimmed": not untrimmed,
        "pool": pool,
    }

    sampleinfo = review_sampleinfo(sampleinfo_file)
    run_creator = Job_Creator(
        log=logger,
        folders=config.folders,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=config.dry,
        config_path=config.config_path,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=sampleinfo,
        run_settings=run_settings,
    )

    ext_refs = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=sampleinfo,
        force=force_update,
    )
    try:
        ext_refs.db_access.check_ref_lock()
    except RefUpdateLockError as e:
        click.echo("ERROR - {}".format(e))
        click.Abort()
    click.echo("INFO - Checking versions of references..")
    try:
        if not skip_update:
            ext_refs.identify_new(project=True)
            ext_refs.update_refs()
            click.echo("INFO - Version check done. Creating sbatch jobs")
        else:
            click.echo("INFO - Skipping version check.")
    except Exception as e:
        click.echo(f"{e}")
    if len(sampleinfo) > 1:
        run_creator.project_job()
    elif len(sampleinfo) == 1:
        run_creator.project_job(single_sample=True)
    else:
        click.Abort()

    done()


@root.group()
@click.pass_context
def utils(ctx):
    """Utilities for specific purposes"""
    pass


@utils.group()
@click.pass_context
def refer(ctx):
    """Manipulates MLST organisms"""
    pass


@utils.command()
@click.argument("sampleinfo_file")
@click.option("--input", help="Full path to project folder", default="")
@click.option(
    "--track",
    help="Run a specific analysis track",
    default="default",
    type=click.Choice(["default", "typing", "qc", "cgmlst"]),
)
@click.option(
    "--dry",
    help="Builds instance without posting to SLURM",
    default=False,
    is_flag=True,
)
@click.option("--email", default="", help="Forced e-mail recipient")
@click.option("--skip_update", default=False, help="Skips downloading of references", is_flag=True)
@click.option(
    "--report",
    default="default",
    type=click.Choice(["default", "typing", "motif_overview", "qc", "json_dump", "st_update"]),
)
@click.option("--output", help="Report output folder", default="")
@pass_config
def finish(
    config: MicroSALTConfig, sampleinfo_file, input, track, dry, email, skip_update, report, output
):
    """Sequence analysis, typing and resistance identification"""
    pool = []
    if email:
        config.regex.mail_recipient = email
    config.dry = dry
    if not os.path.isdir(input):
        click.echo(f"ERROR - Sequence data folder {input} does not exist.")
        click.Abort()
    if output == "":
        output = input
    for subfolder in os.listdir(input):
        if os.path.isdir(f"{input}/{subfolder}"):
            pool.append(subfolder)

    run_settings = {
        "input": input,
        "track": track,
        "dry": dry,
        "email": config.regex.mail_recipient,
        "skip_update": skip_update,
    }

    sampleinfo = review_sampleinfo(sampleinfo_file)
    ext_refs = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=sampleinfo,
    )
    try:
        ext_refs.db_access.check_ref_lock()
    except RefUpdateLockError as e:
        click.echo("ERROR - {}".format(e))
        click.Abort()
    click.echo("INFO - Checking versions of references..")
    try:
        if not skip_update:
            ext_refs.identify_new(project=True)
            ext_refs.update_refs()
            click.echo("INFO - Version check done. Creating sbatch jobs")
        else:
            click.echo("INFO - Skipping version check.")
    except Exception as e:
        click.echo(f"{e}")

    res_scraper = Scraper(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        slurm_header=config.slurm_header,
        regex=config.regex,
        dry=config.dry,
        config_path=config.config_path,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        sampleinfo=sampleinfo,
        input=input,
    )
    if isinstance(sampleinfo, list) and len(sampleinfo) > 1:
        res_scraper.scrape_project()
    else:
        res_scraper.scrape_sample()

    codemonkey = Reporter(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        regex=config.regex,
        sampleinfo=sampleinfo,
        output=output,
        collection=True,
    )
    codemonkey.report(report)
    done()


@refer.command()
@click.argument("organism")
@click.option("--force", help="Redownloads existing organism", default=False, is_flag=True)
@pass_config
def add(config: MicroSALTConfig, organism, force):
    """Adds a new internal organism from pubMLST"""
    referee = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        force=force,
    )
    try:
        referee.add_pubmlst(organism)
    except Exception as e:
        click.echo(e.args[0])
        click.Abort()
    click.echo("INFO - Checking versions of all references..")
    referee = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        force=force,
    )
    referee.update_refs()


@refer.command()
@pass_config
def observe(config: MicroSALTConfig):
    """Lists all stored organisms"""
    refe = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
    )
    click.echo("INFO - Currently stored organisms:")
    for org in sorted(refe.existing_organisms()):
        click.echo(org.replace("_", " ").capitalize())


@utils.command()
@click.argument("sampleinfo_file")
@click.option("--email", default="", help="Forced e-mail recipient")
@click.option(
    "--type",
    default="default",
    type=click.Choice(["default", "typing", "motif_overview", "qc", "json_dump", "st_update"]),
)
@click.option("--output", help="Full path to output folder", default="")
@click.option("--collection", default=False, is_flag=True)
@pass_config
def report(config: MicroSALTConfig, sampleinfo_file, email, type, output, collection):
    """Re-generates report for a project"""
    if email:
        config.regex.mail_recipient = email
    sampleinfo = review_sampleinfo(sampleinfo_file)
    codemonkey = Reporter(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        regex=config.regex,
        sampleinfo=sampleinfo,
        output=output,
        collection=collection,
    )
    codemonkey.report(type)
    done()


@utils.command("get-bigsdb-credentials")
@click.argument("service", type=click.Choice(["pubmlst", "pasteur"]))
@click.option("--species", default=None, help="Species name (required for the 'pasteur' service)")
@pass_config
def get_bigsdb_credentials(config: MicroSALTConfig, service, species):
    """Obtain and store BIGSdb OAuth credentials for SERVICE (pubmlst or pasteur)"""
    get_bigsdb_credentials_main(service, config, species)


@utils.command()
@click.option("--input", help="Full path to project folder", default=os.getcwd())
@click.pass_context
def generate(ctx, input):
    """Creates a blank sample info json for the given input folder"""
    input = os.path.abspath(input)
    project_name = os.path.basename(input)

    defaults = default_sampleinfo.copy()

    pool = []
    if not os.path.isdir(input):
        click.echo(f"ERROR - Sequence data folder {project_name} does not exist.")
        ctx.abort()
    elif input != os.getcwd():
        for subfolder in os.listdir(input):
            if os.path.isdir(f"{input}/{subfolder}"):
                pool.append(defaults.copy())
                pool[-1]["CG_ID_project"] = project_name
                pool[-1]["CG_ID_sample"] = subfolder
    else:
        project_name = "default_sample_info"
        pool.append(defaults.copy())

    with open(f"{os.getcwd()}/{project_name}.json", "w") as output:
        json.dump(pool, output, indent=2)
    click.echo(f"INFO - Created {project_name}.json in current folder")
    done()


@utils.group()
@click.pass_context
def resync(ctx):
    """Updates internal ST with pubMLST equivalent"""


@resync.command()
@click.option(
    "--type",
    default="list",
    type=click.Choice(["report", "list"]),
    help="Output format",
)
@click.option("--customer", default="all", help="Customer id filter")
@click.option("--skip_update", default=False, help="Skips downloading of references", is_flag=True)
@click.option("--email", default="", help="Forced e-mail recipient")
@click.option("--output", help="Full path to output folder", default="")
@pass_config
def review(config: MicroSALTConfig, type, customer, skip_update, email, output):
    """Generates information about novel ST"""
    if email:
        config.regex.mail_recipient = email
    ext_refs = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
    )
    if not skip_update:
        ext_refs.update_refs()
        ext_refs.resync()
    click.echo("INFO - Version check done. Generating output")
    if type == "report":
        codemonkey = Reporter(
            log=logger,
            folders=config.folders,
            threshold=config.threshold,
            regex=config.regex,
            output=output,
        )
        codemonkey.report(type="st_update", customer=customer)
    elif type == "list":
        ext_refs.resync(type=type)
    done()


@resync.command()
@click.option("--force-update", default=False, is_flag=True, help="Forces update")
@pass_config
def update_refs(config: MicroSALTConfig, force_update: bool):
    """Updates all references"""
    ext_refs = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        force=force_update,
    )
    ext_refs.update_refs()
    done()


@resync.command()
@click.option("--force-update", default=False, is_flag=True, help="Forces update")
@pass_config
def update_from_static(config: MicroSALTConfig, force_update: bool):
    """Updates a specific organism"""
    ext_refs = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        force=force_update,
    )
    ext_refs.fetch_external()
    done()


@resync.command()
@click.argument("organism")
@click.option("--force-update", default=False, is_flag=True, help="Forces update")
@click.option("--external", is_flag=True, default=False, help="Updates from external sources")
@pass_config
def update_organism(config: MicroSALTConfig, external: bool, force_update: bool, organism: str):
    """Updates a specific organism"""
    ext_refs = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
        force=force_update,
    )
    ext_refs.update_organism(external=external, organism=organism)
    done()


@resync.command()
@click.argument("sample_name")
@click.option(
    "--force",
    default=False,
    is_flag=True,
    help="Resolves sample without checking for pubMLST match",
)
@pass_config
def overwrite(config: MicroSALTConfig, sample_name, force):
    """Flags sample as resolved"""
    ext_refs = Referencer(
        log=logger,
        folders=config.folders,
        threshold=config.threshold,
        pubmlst=config.pubmlst,
        pasteur=config.pasteur,
        singularity=config.singularity,
        containers=config.containers,
    )
    ext_refs.resync(type="overwrite", sample=sample_name, ignore=force)
    done()
