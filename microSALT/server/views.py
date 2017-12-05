from datetime import date
from flask import Flask, render_template

from microSALT.store.orm_models import Projects, Samples, Seq_types
from microSALT import app

@app.route('/microSALT/')
def start_page():
    projects = Projects.query.all()

    return render_template('start_page.html',
        projects = projects)

@app.route('/microSALT/<project>')
def project_page(project):
    species = list()
    distinct_organisms = Samples.query.filter_by(CG_ID_project=project).distinct()
    for one_guy in distinct_organisms:
        if one_guy.organism not in species:
            species.append(one_guy.organism)
    return render_template('project_page.html',
        organisms = species,
        project = project) 

@app.route('/microSALT/<project>/<organism>')
def report_page(project, organism):
    # Joins all three tables, and displays only hits with 100% hit rate
    sample_info = Samples.query.filter_by(organism=organism).join(Projects).filter_by(CG_ID_project=project).join(Seq_types).all()

    return render_template('report_page.html',
        sample_info = sample_info,
        date = date.today().isoformat())

@app.route('/microSALT/all')
def all_page():
    sample_info = Samples.query.join(Projects).join(Seq_types).all()

    return render_template('debug_page.html',
        sample_info = sample_info,
        date = date.today().isoformat())

