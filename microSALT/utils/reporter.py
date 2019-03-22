"""Generates various reports by tapping into Flask and mySQL
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import requests
import time
import os
import sys
import smtplib

from os.path import basename
from email.mime.text  import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication

from multiprocessing import Process

from microSALT.server.views import app, session, gen_reportdata
from microSALT.store.orm_models import Samples
from microSALT.store.lims_fetcher import LIMS_Fetcher

class Reporter():

  def __init__(self, config, log, name = ""):
    self.name = name
    self.config = config
    self.logger = log
    self.server = Process(target=app.run)
    self.ticketFinder = LIMS_Fetcher(self.config, self.logger)
    self.attachments = list()

  def report(self, type='html'):
    if type == 'html':
      self.gen_html()
    elif type == 'csv':
      self.csv()
    else:
      raise Exception("Report function recieved invalid format")
    self.mail()
    for file in self.attachments:
      os.remove(file)

  def gen_html(self):
    self.start_web()
    name = self.name
    try:
      self.ticketFinder.load_lims_project_info(name)
    except Exception as e:
      self.logger.error("Project {} does not exist".format(name))
      self.kill_flask()
      sys.exit(-1)
    try:
      r = requests.get("http://127.0.0.1:5000/microSALT/{}/all".format(name), allow_redirects=True)
    except Exception as e:
      self.logger.error("Flask instance currently occupied. Possible rogue process. Retry command")
      self.kill_flask()
      sys.exit(-1)
    outname = "{}_microSALT.html".format(self.ticketFinder.data['Customer_ID_project'])
    open(outname, 'wb').write(r.content)
    self.attachments.append(outname)
    self.kill_flask()

  def mail(self):
    """ Requires that the sending PC has a domain"""
    file_name = self.attachments
    msg = MIMEMultipart()
    msg['Subject'] = '{} ({}) Reports'.format(self.name, file_name[0].split('_')[0])
    msg['From'] = 'microSALT'
    msg['To'] = self.config['regex']['mail_recipient']
    
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

  def csv(self):
    name = self.name
    self.ticketFinder.load_lims_project_info(name)
    excel = open("{}.csv".format(name), "w+")
    sample_info = gen_reportdata(name)
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
    inpath = "{}/{}.csv".format(os.getcwd(), name)
    self.attachments.append(inpath)


  def kill_flask(self):
    self.server.terminate()
    self.server.join()
    self.logger.info("Closed webserver on http://127.0.0.1:5000/")
