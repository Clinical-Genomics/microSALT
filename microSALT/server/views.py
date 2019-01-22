from datetime import date
from flask import Flask, render_template
import logging
from io import StringIO, BytesIO

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import *
from sqlalchemy.sql.expression import case, func

from microSALT import config, __version__
from microSALT.store.db_manipulator import app
from microSALT.store.orm_models import Projects, Samples, Seq_types, Versions

engine = create_engine(app.config['SQLALCHEMY_DATABASE_URI'])
Session = sessionmaker(bind=engine)
session = Session()
app.debug = 0
#Removes server start messages
log = logging.getLogger('werkzeug')
log.setLevel(logging.CRITICAL)


@app.route('/')
def start_page():
    projects = session.query(Projects).all()

    return render_template('start_page.html',
        projects = projects)

@app.route('/microSALT/')
def reroute_page():
    projects = session.query(Projects).all()

    return render_template('start_page.html',
        projects = projects)

@app.route('/microSALT/<project>')
def project_page(project):
    organism_groups = list()
    organism_groups.append("all")
    distinct_organisms = session.query(Samples).filter_by(CG_ID_project=project).distinct()
    for one_guy in distinct_organisms:
        if one_guy.organism not in organism_groups and one_guy.organism is not None:
            organism_groups.append(one_guy.organism)
    organism_groups.sort()
    return render_template('project_page.html',
        organisms = organism_groups,
        project = project) 

@app.route('/microSALT/<project>/<organism_group>')
def report_page(project, organism_group):
    sample_info = gen_reportdata(project, organism_group)

    return render_template('report_page.html',
        samples = sample_info['samples'],
        date = date.today().isoformat(),
        version = sample_info['versions'],
        build = __version__)

def gen_reportdata(pid, organism_group='all'):
  """ Queries database for all necessary information for the reports """
  output = dict()
  output['samples'] = list()
  output['versions'] = dict()
  if organism_group=='all':
    sample_info = session.query(Samples).filter(Samples.CG_ID_project==pid)
  else:
    sample_info = session.query(Samples).\
                  filter(Samples.organism==organism_group, Samples.CG_ID_project==project)
  #Sorts sample names
  sample_info = sorted(sample_info, key=lambda sample: \
                int(sample.CG_ID_sample.replace(sample.CG_ID_project, '')[1:]))
  for s in sample_info:
    s.ST_status=str(s.ST)
    if s.Customer_ID_sample.startswith('NTC') or s.Customer_ID_sample.startswith('0-') or s.Customer_ID_sample.startswith('NK-') or s.Customer_ID_sample.startswith('NEG') or s.Customer_ID_sample.startswith('CTRL') or s.Customer_ID_sample.startswith('Neg'):
      s.ST_status = 'Control (prefix)'
    elif s.ST < 0:
      if s.ST == -1:
        s.ST_status = 'Unavailable'
      elif s.ST == -4:
        s.ST_status = 'Novel'
      else:
        s.ST_status='None'

    if 'Control' in s.ST_status:
      s.threshold = 'Passed'
    elif s.ST == -3:
      s.threshold = 'Failed'
    elif hasattr(s, 'seq_types') and s.seq_types != [] or s.ST == -2:
      near_hits=0
      s.threshold = 'Passed'
      for seq_type in s.seq_types:
        #Identify single deviating allele
        if seq_type.st_predictor and seq_type.identity >= config["threshold"]["mlst_novel_id"] and config["threshold"]["mlst_id"] > seq_type.identity and 1-abs(1-seq_type.span) >= config["threshold"]["mlst_span"]:
          near_hits = near_hits + 1
        elif (seq_type.identity < config["threshold"]["mlst_novel_id"] or seq_type.span < config["threshold"]["mlst_span"]) and seq_type.st_predictor:
          s.threshold = 'Failed'

      if near_hits > 0 and s.threshold == 'Passed':
        s.ST_status = 'Novel ({} alleles)'.format(near_hits)
    else:
      s.threshold = 'Failed'

    #Resistence filter
    for r in s.resistances:
      if (s.ST > 0 or 'Novel' in s.ST_status ) and (r.identity >= config["threshold"]["res_id"] and r.span >= config["threshold"]["res_span"]) or (s.ST < 0 and not 'Novel' in s.ST_status):
        r.threshold = 'Passed'
      else:
        r.threshold = 'Failed'

    #Seq_type and resistance sorting
    s.seq_types=sorted(s.seq_types, key=lambda x: x.loci)
    s.resistances=sorted(s.resistances, key=lambda x: x.instance)
    output['samples'].append(s)

  versions = session.query(Versions).all()
  for version in versions:
    name = version.name[8:]
    output['versions'][name] = version.version

  return output
