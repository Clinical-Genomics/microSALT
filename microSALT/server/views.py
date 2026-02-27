import logging
import subprocess
from datetime import date
from pathlib import Path

from jinja2 import Environment, FileSystemLoader

from microSALT import __version__, preset_config
from microSALT.store.orm_models import (
    Collections,
    Reports,
    Samples,
    Versions,
)

from microSALT.store.database import get_session

# Removes server start messages
log = logging.getLogger("werkzeug")
log.setLevel(logging.CRITICAL)

TEMPLATE_FOLDER = Path(__file__).parent / "templates"


def _make_url_for(project=None):
    """Creates a url_for stub for use in static Jinja2 templates."""

    def url_for(endpoint, **kwargs):
        if endpoint == "static":
            return kwargs.get("filename", "")
        elif endpoint == "start_page":
            return "#"
        elif endpoint == "project_page":
            p = kwargs.get("project", project or "")
            return "#project-{}".format(p)
        elif endpoint == "typing_page":
            p = kwargs.get("project", project or "")
            og = kwargs.get("organism_group", "all")
            return "#typing-{}-{}".format(p, og)
        return "#"

    return url_for


def render_template(template_folder, template_name, **context):
    """Renders a template using Jinja2 directly to avoid Flask overhead"""
    template_loader = FileSystemLoader(searchpath=str(template_folder))
    jinja_env = Environment(loader=template_loader)
    template = jinja_env.get_template(template_name)
    if "url_for" not in context:
        context["url_for"] = _make_url_for()
    return template.render(**context)


def project_page(project, template_folder: Path = TEMPLATE_FOLDER):
    session = get_session()
    organism_groups = list()
    organism_groups.append("all")
    distinct_organisms = session.query(Samples).filter_by(CG_ID_project=project).distinct()
    for one_guy in distinct_organisms:
        if one_guy.organism not in organism_groups and one_guy.organism is not None:
            organism_groups.append(one_guy.organism)
    organism_groups.sort()
    return render_template(
        template_folder=template_folder,
        template_name="project_page.html",
        organisms=organism_groups,
        project=project,
        url_for=_make_url_for(project),
    )


def alignment_page(project, template_folder: Path = TEMPLATE_FOLDER):
    sample_info = gen_reportdata(project)

    return render_template(
        template_folder=template_folder,
        template_name="alignment_page.html",
        samples=sample_info["samples"],
        topsample=sample_info["single_sample"],
        date=date.today().isoformat(),
        version=sample_info["versions"],
        user=sample_info["user"],
        threshold=preset_config["threshold"],
        reports=sample_info["reports"],
        build=__version__,
    )


def render_alignment_page(project, template_folder: Path = TEMPLATE_FOLDER):
    return alignment_page(project, template_folder=template_folder)


def typing_page(project, organism_group, template_folder: Path = TEMPLATE_FOLDER):
    sample_info = gen_reportdata(project, organism_group)

    return render_template(
        template_folder=template_folder,
        template_name="typing_page.html",
        samples=sample_info["samples"],
        topsample=sample_info["single_sample"],
        date=date.today().isoformat(),
        version=sample_info["versions"],
        user=sample_info["user"],
        threshold=preset_config["threshold"],
        verified_organisms=preset_config["regex"]["verified_organisms"],
        reports=sample_info["reports"],
        build=__version__,
    )


def render_typing_page(project, organism_group, template_folder: Path = TEMPLATE_FOLDER):
    return typing_page(project, organism_group, template_folder=template_folder)


def STtracker_page(customer, template_folder: Path = TEMPLATE_FOLDER):
    sample_info = gen_reportdata(project_id="all", organism_group="all")
    final_samples = list()
    for s in sample_info["samples"]:
        if customer == "all" or s.projects.Customer_ID == customer:
            if s.pubmlst_ST != -1 and s.ST < 0:
                final_samples.append(s)

    final_samples = sorted(final_samples, key=lambda sample: (sample.CG_ID_sample))

    return render_template(
        template_folder=template_folder,
        template_name="STtracker_page.html",
        date=date.today().isoformat(),
        internal=final_samples,
    )


def gen_collectiondata(collect_id=[]):
    """Queries database using a set of samples"""
    session = get_session()
    samples = session.query(Collections).filter(Collections.ID_collection == collect_id).all()
    sample_ids = [s.CG_ID_sample for s in samples]
    sample_info = session.query(Samples).filter(Samples.CG_ID_sample.in_(sample_ids))
    sample_info = gen_add_info(sample_info)
    return sample_info


def gen_reportdata(project_id="all", organism_group="all"):
    """Queries database for all necessary information for the reports"""
    session = get_session()
    if project_id == "all" and organism_group == "all":
        sample_info = session.query(Samples)
    elif project_id == "all":
        sample_info = session.query(Samples).filter(Samples.organism == organism_group)
    elif organism_group == "all":
        sample_info = session.query(Samples).filter(Samples.CG_ID_project == project_id)
    else:
        sample_info = session.query(Samples).filter(
            Samples.CG_ID_project == project_id, Samples.organism == organism_group
        )

    sample_info = gen_add_info(sample_info)

    reports = session.query(Reports).filter(Reports.CG_ID_project == project_id).all()
    sample_info["reports"] = reports = sorted(reports, key=lambda x: x.version, reverse=True)

    return sample_info


def gen_add_info(sample_info=dict()):
    """Enhances a sample info struct by adding ST_status, threshold info, versioning and sorting"""
    session = get_session()
    # Set ST status
    output = dict()
    output["samples"] = list()
    output["versions"] = dict()
    output["single_sample"] = ""

    # Sorts sample names
    valid = True
    for sam in sample_info.all():
        if sam.CG_ID_project is None:
            valid = False
            break
    if valid:
        try:
            sample_info = sorted(
                sample_info,
                key=lambda sample: int(sample.CG_ID_sample.replace(sample.CG_ID_project, "")[1:]),
            )
        except ValueError:
            pass

    for s in sample_info:
        s.CG_ID_project = s.projects.CG_ID_project
        s.ST_status = str(s.ST)
        if s.Customer_ID_sample is not None:
            if (
                s.Customer_ID_sample.startswith("NTC")
                or s.Customer_ID_sample.startswith("0-")
                or s.Customer_ID_sample.startswith("NK-")
                or s.Customer_ID_sample.startswith("NEG")
                or s.Customer_ID_sample.startswith("CTRL")
                or s.Customer_ID_sample.startswith("Neg")
                or s.Customer_ID_sample.startswith("blank")
                or s.Customer_ID_sample.startswith("dual-NTC")
            ):
                s.ST_status = "Kontroll (prefix)"

        if "Kontroll" in s.ST_status or "Control" in s.ST_status or s.ST == -1:
            s.threshold = "-"
        elif s.ST == -3:
            s.threshold = "Failed"
        elif hasattr(s, "seq_types") and s.seq_types != [] or s.ST == -2:
            near_hits = 0
            s.threshold = "Passed"
            for seq_type in s.seq_types:
                # Identify single deviating allele
                if (
                    seq_type.st_predictor
                    and seq_type.identity >= preset_config["threshold"]["mlst_novel_id"]
                    and preset_config["threshold"]["mlst_id"] > seq_type.identity
                    and 1 - abs(1 - seq_type.span)
                    >= (preset_config["threshold"]["mlst_span"] / 100.0)
                ):
                    near_hits = near_hits + 1
                elif (
                    seq_type.identity < preset_config["threshold"]["mlst_novel_id"]
                    or seq_type.span < (preset_config["threshold"]["mlst_span"] / 100.0)
                ) and seq_type.st_predictor:
                    s.threshold = "Failed"

            if near_hits > 0 and s.threshold == "Passed":
                s.ST_status = f"Okänd ({near_hits} allele[r])"
        else:
            s.threshold = "Failed"

        if not ("Control" in s.ST_status or "Kontroll" in s.ST_status) and s.ST < 0:
            if s.ST == -1:
                s.ST_status = "Data saknas"
            elif s.ST <= -4 or s.ST == -2:
                s.ST_status = "Okänd (Novel ST, Novel allele[r])"
            else:
                s.ST_status = "None"

        # Resistence filter
        for r in s.resistances:
            if (
                r.identity >= preset_config["threshold"]["motif_id"]
                and r.span >= preset_config["threshold"]["motif_span"] / 100.0
            ):
                r.threshold = "Passed"
            else:
                r.threshold = "Failed"
        for v in s.expacs:
            if (
                v.identity >= preset_config["threshold"]["motif_id"]
                and v.span >= preset_config["threshold"]["motif_span"] / 100.0
            ):
                v.threshold = "Passed"
            else:
                v.threshold = "Failed"

        # Seq_type and resistance sorting
        s.seq_types = sorted(s.seq_types, key=lambda x: x.loci)
        s.resistances = sorted(s.resistances, key=lambda x: x.instance)
        s.expacs = sorted(s.expacs, key=lambda x: x.gene)
        output["samples"].append(s)
        output["single_sample"] = s

    versions = session.query(Versions).all()
    for version in versions:
        name = version.name[8:]
        output["versions"][name] = version.version

    process = subprocess.Popen("id -un".split(), stdout=subprocess.PIPE)
    user, error = process.communicate()
    output["user"] = user.decode("utf-8").replace(".", " ").title()

    return output
