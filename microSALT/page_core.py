import os
import yaml
from flask import Flask, render_template
from flask_sqlalchemy import SQLAlchemy
from microSALT import app, db

from microSALT.tables.samples import Samples
from microSALT.tables.seq_types import Seq_types
from microSALT.tables.profiles import Profiles
from microSALT.db_manipulator import DB_Manipulator

@app.route('/microSALT/')
def start_page():
    samples = Samples.query.all()
    return render_template('start_page.html',
        samples = samples)


@app.route('/microSALT/<sample>')
def sample_page(sample):
    seq_types = Seq_types.query.filter_by(CG_ID_sample=sample)
    samples = Samples.query.filter_by(CG_ID_sample=sample).first()
    reduced_seqtypes = Seq_types.query.filter(Seq_types.CG_ID_sample==sample , Seq_types.identity==100)


    return render_template('sample_page.html',
        seq_types = reduced_seqtypes,
        samples = samples) 
