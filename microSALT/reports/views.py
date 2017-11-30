from flask import Flask, render_template

from microSALT.tables.samples import Samples
from microSALT.tables.seq_types import Seq_types
from microSALT import app

from datetime import date

@app.route('/microSALT/')
def start_page():
    samples = Samples.query.all()
    projects = list(set([(samp.CG_ID_project,samp.Customer_ID_project) for samp in samples]))

    return render_template('start_page.html',
        projects = projects)

@app.route('/microSALT/<project>')
def project_page(project):
    organisms = []
    #TODO: Establish organism by database query
    all_organisms = ['enterococcus_faecalis','enterococcus_faecium','escherichia_coli','klebsiella_pneumoniae','staphylococcus_aureus']
    for organism in all_organisms:
        samples = Samples.query.filter_by(organism=organism, CG_ID_project=project).all()
        if samples:
            organisms.append(organism)

    return render_template('project_page.html',
        organisms = organisms,
        project = project) 

@app.route('/microSALT/<project>/<organism>')
def report_page(project, organism):
    samples = Samples.query.filter_by(organism=organism, CG_ID_project=project).all()
    reduced_samples = [ {'sample':sample, 'seq_types' : Seq_types.query.filter_by(CG_ID_sample=sample.CG_ID_sample, identity=100)} for sample in samples]

    print date.today().isoformat()
    return render_template('report_page.html',
        reduced_samples = reduced_samples,
        date = date.today().isoformat())
