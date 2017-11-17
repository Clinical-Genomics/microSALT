import os
import yaml
from flask import Flask, render_template
from flask_sqlalchemy import SQLAlchemy
from microSALT import app, db

class Blast(db.Model):
    __tablename__ = 'blast'
    run = db.Column(db.String(40), primary_key=True)
    date_analysis = db.Column(db.DateTime)
    organism = db.Column(db.String(30))
    loci = db.Column(db.String(10))
    assumed_ST = db.Column(db.SmallInteger)
    allele = db.Column(db.SmallInteger)
    haplotype = db.Column(db.String(5))
    contig_name = db.Column(db.String(20), primary_key=True)
    contig_length = db.Column(db.Integer)
    contig_coverage = db.Column( db.Float(6,2))
    identity = db.Column(db.Float(3,2))
    evalue = db.Column(db.String(10))
    bitscore = db.Column(db.SmallInteger)
    contig_start = db.Column(db.Integer)
    contig_end = db.Column(db.Integer)
    loci_start = db.Column(db.Integer)
    loci_end = db.Column(db.Integer)


@app.route('/microSALT/')
def start_page():
    blast = Blast.query.all()
    return render_template('start_page.html',
        blast  = blast)


@app.route('/microSALT/<sample>')
def sample_page(sample):
    blast = Blast.query.filter_by(run=sample)
    return render_template('sample_page.html',
        blast  = blast,
        sample_name = sample,
        example_sample = blast[0])
