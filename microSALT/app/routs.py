import os
import sys
from flask import Flask, render_template
from flask_sqlalchemy import SQLAlchemy
from microSALT import db
from microSALT.tables.samples import Samples
#from microSALT.app import app

@app.route('/microSALT/')
def start_page():
    blast = Samples.query.all()
    return render_template('start_page.html',
        blast  = blast)


@app.route('/microSALT/<sample>')
def sample_page(sample):
    blast = Samples.query.filter_by(CG_ID_sample=sample)
    #print dir(blast[0])
    return render_template('sample_page.html',
        blast  = blast,
        sample_name = sample,
        example_sample = blast[0])
