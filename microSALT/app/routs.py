import os
import sys
from flask import Flask, render_template
from flask_sqlalchemy import SQLAlchemy
from bla bla import app
from blo blo import db, Blast som ska ersättas med alla modellerna

@app.route('/microSALT/')
def start_page():
    blast = Blast.query.all()
    return render_template('start_page.html',
        blast  = blast)


@app.route('/microSALT/<sample>')
def sample_page(sample):
    blast = Blast.query.filter_by(run=sample)
    print dir(blast[0])
    return render_template('sample_page.html',
        blast  = blast,
        sample_name = sample,
        example_sample = blast[0])
