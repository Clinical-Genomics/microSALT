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
from microSALT.store.database import get_scoped_session_registry, initialize_database
from microSALT.utils.job_creator import Job_Creator
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

    log_dir = os.path.dirname(config.folders.log_file)
    if log_dir and not pathlib.Path(log_dir).exists():
        os.makedirs(log_dir)
    proc = subprocess.Popen(
        f"touch {config.folders.log_file}".split(), stdout=subprocess.PIPE
    )
    proc.communicate()

    folder_paths = [
        config.folders.results,
        config.folders.reports,
        config.folders.seqdata,
        config.folders.profiles,
        config.folders.references,
        config.folders.resistances,
        config.folders.genomes,
        config.folders.credentials,
        config.folders.adapters,
    ]
    for path in folder_paths:
        p = pathlib.Path(os.path.expandvars(os.path.expanduser(path)))
        if not p.exists():
            os.makedirs(p)


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
    _ensure_directories(cfg)
    initialize_database(cfg.database.SQLALCHEMY_DATABASE_URI)
    setup_logger(logging_level=logging_level, log_file=cfg.folders.log_file)
    logger.setLevel(logging_levels[logging_level])
    for handler in logger.handlers:
        handler.setLevel(logging_levels[logging_level])
    ctx.obj = {}
    ctx.obj["config"] = cfg
    ctx.call_on_close(teardown_session)


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
@click.pass_context
def analyse(
    ctx,
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
        ctx.obj["config"].regex.mail_recipient = email
    ctx.obj["config"].dry = dry
    if not os.path.isdir(input):
        click.echo(f"ERROR - Sequence data folder {input} does not exist.")
        ctx.abort()
    for subfolder in os.listdir(input):
        if os.path.isdir(f"{input}/{subfolder}"):
            pool.append(subfolder)

    run_settings = {
        "input": input,
        "dry": dry,
        "email": ctx.obj["config"].regex.mail_recipient,
        "skip_update": skip_update,
        "trimmed": not untrimmed,
        "pool": pool,
    }

    sampleinfo = review_sampleinfo(sampleinfo_file)
    cfg = ctx.obj["config"]
    run_creator = Job_Creator(
        log=logger,
        folders=cfg.folders,
        slurm_header=cfg.slurm_header,
        regex=cfg.regex,
        dry=cfg.dry,
        config_path=cfg.config_path,
        threshold=cfg.threshold,
        pubmlst=cfg.pubmlst,
        pasteur=cfg.pasteur,
        singularity=cfg.singularity,
        containers=cfg.containers,
        sampleinfo=sampleinfo,
        run_settings=run_settings,
    )

    ext_refs = Referencer(
        log=logger,
        folders=cfg.folders,
        threshold=cfg.threshold,
        pubmlst=cfg.pubmlst,
        pasteur=cfg.pasteur,
        singularity=cfg.singularity,
        containers=cfg.containers,
        sampleinfo=sampleinfo,
        force=force_update,
    )
    try:
        ext_refs.db_access.check_ref_lock()
    except RefUpdateLockError as e:
        click.echo("ERROR - {}".format(e))
        ctx.abort()
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
        ctx.abort()

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
@click.pass_context
def finish(ctx, sampleinfo_file, input, track, dry, email, skip_update, report, output):
    """Sequence analysis, typing and resistance identification"""
    pool = []
    if email:
        ctx.obj["config"].regex.mail_recipient = email
    ctx.obj["config"].dry = dry
    if not os.path.isdir(input):
        click.echo(f"ERROR - Sequence data folder {input} does not exist.")
        ctx.abort()
    if output == "":
        output = input
    for subfolder in os.listdir(input):
        if os.path.isdir(f"{input}/{subfolder}"):
            pool.append(subfolder)

    run_settings = {
        "input": input,
        "track": track,
        "dry": dry,
        "email": ctx.obj["config"].regex.mail_recipient,
        "skip_update": skip_update,
    }

    sampleinfo = review_sampleinfo(sampleinfo_file)
    cfg = ctx.obj["config"]
    ext_refs = Referencer(
        log=logger,
        folders=cfg.folders,
        threshold=cfg.threshold,
        pubmlst=cfg.pubmlst,
        pasteur=cfg.pasteur,
        singularity=cfg.singularity,
        containers=cfg.containers,
        sampleinfo=sampleinfo,
    )
    try:
        ext_refs.db_access.check_ref_lock()
    except RefUpdateLockError as e:
        click.echo("ERROR - {}".format(e))
        ctx.abort()
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
        folders=cfg.folders,
        threshold=cfg.threshold,
        slurm_header=cfg.slurm_header,
        regex=cfg.regex,
        dry=cfg.dry,
        config_path=cfg.config_path,
        pubmlst=cfg.pubmlst,
        pasteur=cfg.pasteur,
        singularity=cfg.singularity,
        containers=cfg.containers,
        sampleinfo=sampleinfo,
        input=input,
    )
    if isinstance(sampleinfo, list) and len(sampleinfo) > 1:
        res_scraper.scrape_project()
    else:
        res_scraper.scrape_sample()

    codemonkey = Reporter(
        log=logger,
        folders=cfg.folders,
        threshold=cfg.threshold,
        regex=cfg.regex,
        sampleinfo=sampleinfo,
        output=output,
        collection=True,
    )
    codemonkey.report(report)
    done()


@refer.command()
@click.argument("organism")
@click.option("--force", help="Redownloads existing organism", default=False, is_flag=True)
@click.pass_context
def add(ctx, organism, force):
    """Adds a new internal organism from pubMLST"""
    cfg = ctx.obj["config"]
    referee = Referencer(
        log=logger, folders=cfg.folders, threshold=cfg.threshold,
        pubmlst=cfg.pubmlst, pasteur=cfg.pasteur,
        singularity=cfg.singularity, containers=cfg.containers, force=force,
    )
    try:
        referee.add_pubmlst(organism)
    except Exception as e:
        click.echo(e.args[0])
        ctx.abort()
    click.echo("INFO - Checking versions of all references..")
    referee = Referencer(
        log=logger, folders=cfg.folders, threshold=cfg.threshold,
        pubmlst=cfg.pubmlst, pasteur=cfg.pasteur,
        singularity=cfg.singularity, containers=cfg.containers, force=force,
    )
    referee.update_refs()


@refer.command()
@click.pass_context
def observe(ctx):
    """Lists all stored organisms"""
    cfg = ctx.obj["config"]
    refe = Referencer(
        log=logger, folders=cfg.folders, threshold=cfg.threshold,
        pubmlst=cfg.pubmlst, pasteur=cfg.pasteur,
        singularity=cfg.singularity, containers=cfg.containers,
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
@click.pass_context
def report(ctx, sampleinfo_file, email, type, output, collection):
    """Re-generates report for a project"""
    if email:
        ctx.obj["config"].regex.mail_recipient = email
    cfg = ctx.obj["config"]
    sampleinfo = review_sampleinfo(sampleinfo_file)
    codemonkey = Reporter(
        log=logger,
        folders=cfg.folders,
        threshold=cfg.threshold,
        regex=cfg.regex,
        sampleinfo=sampleinfo,
        output=output,
        collection=collection,
    )
    codemonkey.report(type)
    done()


@utils.command()
@click.pass_context
def view(ctx):
    """Starts an interactive webserver for viewing"""
    cfg = ctx.obj["config"]
    codemonkey = Reporter(log=logger, folders=cfg.folders, threshold=cfg.threshold, regex=cfg.regex)
    codemonkey.start_web()


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
@click.pass_context
def review(ctx, type, customer, skip_update, email, output):
    """Generates information about novel ST"""
    if email:
        ctx.obj["config"].regex.mail_recipient = email
    cfg = ctx.obj["config"]
    ext_refs = Referencer(
        log=logger, folders=cfg.folders, threshold=cfg.threshold,
        pubmlst=cfg.pubmlst, pasteur=cfg.pasteur,
        singularity=cfg.singularity, containers=cfg.containers,
    )
    if not skip_update:
        ext_refs.update_refs()
        ext_refs.resync()
    click.echo("INFO - Version check done. Generating output")
    if type == "report":
        codemonkey = Reporter(log=logger, folders=cfg.folders, threshold=cfg.threshold, regex=cfg.regex, output=output)
        codemonkey.report(type="st_update", customer=customer)
    elif type == "list":
        ext_refs.resync(type=type)
    done()


@resync.command()
@click.option("--force-update", default=False, is_flag=True, help="Forces update")
@click.pass_context
def update_refs(ctx, force_update: bool):
    """Updates all references"""
    cfg = ctx.obj["config"]
    ext_refs = Referencer(
        log=logger, folders=cfg.folders, threshold=cfg.threshold,
        pubmlst=cfg.pubmlst, pasteur=cfg.pasteur,
        singularity=cfg.singularity, containers=cfg.containers, force=force_update,
    )
    ext_refs.update_refs()
    done()


@resync.command()
@click.option("--force-update", default=False, is_flag=True, help="Forces update")
@click.pass_context
def update_from_static(ctx, force_update: bool):
    """Updates a specific organism"""
    cfg = ctx.obj["config"]
    ext_refs = Referencer(
        log=logger, folders=cfg.folders, threshold=cfg.threshold,
        pubmlst=cfg.pubmlst, pasteur=cfg.pasteur,
        singularity=cfg.singularity, containers=cfg.containers, force=force_update,
    )
    ext_refs.fetch_external()
    done()


@resync.command()
@click.argument("organism")
@click.option("--force-update", default=False, is_flag=True, help="Forces update")
@click.option("--external", is_flag=True, default=False, help="Updates from external sources")
@click.pass_context
def update_organism(ctx, external: bool, force_update: bool, organism: str):
    """Updates a specific organism"""
    cfg = ctx.obj["config"]
    ext_refs = Referencer(
        log=logger, folders=cfg.folders, threshold=cfg.threshold,
        pubmlst=cfg.pubmlst, pasteur=cfg.pasteur,
        singularity=cfg.singularity, containers=cfg.containers, force=force_update,
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
@click.pass_context
def overwrite(ctx, sample_name, force):
    """Flags sample as resolved"""
    cfg = ctx.obj["config"]
    ext_refs = Referencer(
        log=logger, folders=cfg.folders, threshold=cfg.threshold,
        pubmlst=cfg.pubmlst, pasteur=cfg.pasteur,
        singularity=cfg.singularity, containers=cfg.containers,
    )
    ext_refs.resync(type="overwrite", sample=sample_name, ignore=force)
    done()
