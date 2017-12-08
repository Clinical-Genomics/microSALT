from datetime import date
from flask import Flask, render_template

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker, aliased
from microSALT.store.orm_models import Projects, Samples, Seq_types
from microSALT import app

engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'])
Session = sessionmaker(bind=engine)
session = Session()

@app.route('/microSALT/')
def start_page():
    projects = session.query(Projects).all()

    return render_template('start_page.html',
        projects = projects)

@app.route('/microSALT/<project>')
def project_page(project):
    species = list()
    distinct_organisms = session.query(Samples).filter_by(CG_ID_project=project).distinct()
    for one_guy in distinct_organisms:
        if one_guy.organism not in species:
            species.append(one_guy.organism)
    return render_template('project_page.html',
        organisms = species,
        project = project) 

@app.route('/microSALT/<project>/<organism>')
def report_page(project, organism):
    #Wild wacky foreign key adventures
    sample_info = session.query(Samples).filter(Samples.organism==organism, Samples.CG_ID_project==project).all()

    return render_template('report_page.html',
        samples = sample_info,
        date = date.today().isoformat())
