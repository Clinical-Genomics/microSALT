"""Generates various reports by tapping into Flask and mySQL
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import json
import requests
import time
import os
import sys
import smtplib

from datetime import datetime

from os.path import basename
from email.mime.text  import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication

from multiprocessing import Process

from microSALT.server.views import app, session, gen_reportdata
from microSALT.store.orm_models import Samples
from microSALT.store.lims_fetcher import LIMS_Fetcher

class Reporter():

  def __init__(self, config, log, name = "", output = ""):
    self.name = name
    if output == "":
      self.output = os.getcwd()
    else:
      self.output = output
    self.config = config
    self.logger = log
    for k, v in config.items():
      app.config[k] = v
    self.server = Process(target=app.run)
    self.ticketFinder = LIMS_Fetcher(self.config, self.logger)
    self.attachments = list()
    self.error = False
    self.dt = datetime.now() 
    self.now = time.strftime("{}.{}.{}_{}.{}.{}".\
    format(self.dt.year, self.dt.month, self.dt.day, self.dt.hour, self.dt.minute, self.dt.second))

  def report(self, type='default', customer='all'):
    self.start_web()
    if type == 'json_dump':
      self.gen_json()
    elif type == 'default':
      self.gen_typing()
      self.gen_qc()
      #self.gen_resistence()
    elif type == 'typing':
      self.gen_typing()
    elif type == 'resistance_overview':
      self.gen_resistence()
    elif type == 'qc':
      self.gen_qc()
    elif type == 'st_update':
      self.gen_STtracker(customer)
    else:
      raise Exception("Report function recieved invalid format")
    self.kill_flask()
    self.mail()
    if self.output == "" or self.output == os.getcwd():
      for file in self.attachments:
        os.remove(file)

  def gen_STtracker(self, customer="all"):
    self.name ="Sequence Type Update"
    try:
      r = requests.get("http://127.0.0.1:5000/microSALT/STtracker/{}".format(customer), allow_redirects=True)
    except Exception as e:
      self.logger.error("Flask instance currently occupied. Possible rogue process. Retry command")
      self.kill_flask()
      sys.exit(-1)
    outname = "{}/ST_updates_{}.html".format(self.output, self.now)
    open(outname, 'wb').write(r.content)
    self.attachments.append(outname)

  def gen_qc(self):
    try:
      self.ticketFinder.load_lims_project_info(self.name)
    except Exception as e:
      self.logger.error("Project {} does not exist".format(self.name))
      self.kill_flask()
      sys.exit(-1)
    try:
      q = requests.get("http://127.0.0.1:5000/microSALT/{}/qc".format(self.name), allow_redirects=True)
      outqc = "{}/{}_QC_{}.html".format(self.output, self.ticketFinder.data['Customer_ID_project'], self.now)
      open(outqc, 'wb').write(q.content)
      self.attachments.append(outqc)  
    except Exception as e:
      self.logger.error("Flask instance currently occupied. Possible rogue process. Retry command")
      self.error = True
 
  def gen_typing(self):
    try:
      self.ticketFinder.load_lims_project_info(self.name)
    except Exception as e:
      self.logger.error("Project {} does not exist".format(self.name))
      self.kill_flask()
      sys.exit(-1)
    try:
      r = requests.get("http://127.0.0.1:5000/microSALT/{}/typing/all".format(self.name), allow_redirects=True)
      outtype = "{}/{}_Typing_{}.html".format(self.output, self.ticketFinder.data['Customer_ID_project'], self.now)
      open(outtype, 'wb').write(r.content)
      self.attachments.append(outtype)
    except Exception as e:
      self.logger.error("Flask instance currently occupied. Possible rogue process. Retry command")
      self.error = True

  def gen_resistence(self):
    self.ticketFinder.load_lims_project_info(self.name)
    output = "{}/{}_{}.csv".format(self.output,self.name,self.now)
    excel = open(output, "w+")
    sample_info = gen_reportdata(self.name)
    resdict = dict()

    #Load ALL resistance & genes
    for s in sample_info['samples']:
      for r in s.resistances:
        if not (r.resistance in resdict.keys()) and r.threshold == 'Passed':
          resdict[r.resistance] =list()
        if r.threshold == 'Passed' and not r.gene in resdict[r.resistance]:
          resdict[r.resistance].append(r.gene)
    for k, v in resdict.items():
      resdict[k] = sorted(v)

    #Header
    topline = ",,,"
    botline = "CG Sample ID,Sample ID,Organism,Sequence Type,Thresholds"
    for k in sorted(resdict.keys()):
      genes = [''] * len(resdict[k])
      topline += ",,{}{}".format(k.replace(',',';'),','.join(genes))
      botline += ",,{}".format(','.join(sorted(resdict[k])))
    excel.write("{}\n".format(topline))
    excel.write("{}\n".format(botline))

    #Individual searches
    for s in sample_info['samples']:
      tdict = dict()
      pref = "{},{},{},{},{}".format(s.CG_ID_sample,s.Customer_ID_sample, s.organism, s.ST_status, s.threshold)
      #Load single sample
      for r in s.resistances:
        if not (r.resistance in tdict.keys()) and r.threshold == 'Passed':
          tdict[r.resistance] =list()
        if r.threshold == 'Passed' and not r.gene in tdict[r.resistance]:
          tdict[r.resistance].append(r.gene)
      #Compare single sample to all
      hits = ""
      for res in sorted(resdict.keys()):
        if res in tdict.keys():
          hits += ",1"
          for gen in sorted(resdict[res]):
            hits += ","
            if gen in tdict[res]:
              hits +="1"
        else:
          #Commas eq to res + gen length
          hits += ",,"
          pad = [''] * len(resdict[res])
          hits += ','.join(pad)

      excel.write("{}{}\n".format(pref, hits))

    excel.close()
    self.attachments.append(output)

  def gen_json(self):
    report = dict()
    output = "{}/{}_{}.json".format(self.output, self.name, self.now)

    sample_info = gen_reportdata(self.name)
    analyses = ['blast_pubmlst', 'quast_assembly', 'blast_resfinder_resistence', 'picard_markduplicate', 'microsalt_samtools_stats']
    for s in sample_info['samples']:
      report[s.CG_ID_sample] = dict()
      for a in analyses:
        if a == 'blast_resfinder_resistence':
          report[s.CG_ID_sample][a] = list()
        else:
          report[s.CG_ID_sample][a] = dict()

      report[s.CG_ID_sample]['blast_pubmlst'] = {'sequence_type':s.ST_status, 'thresholds':s.threshold}
      report[s.CG_ID_sample]['quast_assembly'] = {'estimated_genome_length':s.genome_length, 'gc_percentage':float(s.gc_percentage), 'n50':s.n50, 'necessary_contigs':s.contigs}

      for r in s.resistances:
        if not (r.gene in report[s.CG_ID_sample]['blast_resfinder_resistence']) and r.threshold == 'Passed':
          report[s.CG_ID_sample]['blast_resfinder_resistence'].append(r.gene)

    #json.dumps(report) #Dumps the json directly
    with open(output, 'w') as outfile:
      json.dump(report, outfile)
    self.attachments.append(output)

  def mail(self):
    file_name = self.attachments
    msg = MIMEMultipart()
    if not self.error:
      msg['Subject'] = '{} ({}) Reports'.format(self.name, file_name[0].split('_')[0])
    else:
      msg['Subject'] = '{} ({}) Failed Generating Report'.format(self.name, file_name[0].split('_')[0])
    msg['From'] = 'microSALT'
    msg['To'] = self.config['regex']['mail_recipient']
   
    if not self.error: 
      for file in self.attachments:
        part = MIMEApplication(open(file).read())  
        part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(file))
        msg.attach(part)

    s = smtplib.SMTP('localhost')
    s.connect()
    s.sendmail(msg['From'], msg['To'], msg.as_string())
    s.quit()
    self.logger.info("Mail containing report sent to {}".format(msg['To'])) 

  def start_web(self):
    self.server.start()
    self.logger.info("Started webserver on http://127.0.0.1:5000/")
    #Hinders requests before server goes up
    time.sleep(0.05)

  def kill_flask(self):
    self.server.terminate()
    self.server.join()
    self.logger.info("Closed webserver on http://127.0.0.1:5000/")
