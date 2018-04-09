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

  def report(self):
    self.gen_html()
    self.gen_csv()
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
    r = requests.get("http://127.0.0.1:5000/microSALT/{}/all".format(name), allow_redirects=True)
    outname = "{}_microSALT.html".format(self.ticketFinder.data['Customer_ID_project'])
    open(outname, 'wb').write(r.content)
    self.kill_flask()
    self.attachments.append(outname)

  def mail(self):
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

  def gen_csv(self):
    name = self.name
    self.ticketFinder.load_lims_project_info(name)
    excel = open("{}.csv".format(name), "w+")
    sample_info = gen_reportdata(name)
    excel.write("Customer_ID_sample,CG_ID_sample,organism,ST,Thresholds\n")
    for s in sample_info['samples']:
      if s.organism != None:
        organism_fix = s.organism.replace('_', ' ').capitalize()
      else:
        organism_fix= "None"
      excel.write("{},{},{},{},{}\n".format(s.Customer_ID_sample, s.CG_ID_sample,\
                    organism_fix, s.ST_status,s.threshold))
    excel.close()
    inpath = "{}/{}.csv".format(os.getcwd(), name)
    self.attachments.append(inpath)

  def kill_flask(self):
    self.server.terminate()
    self.server.join()
    self.logger.info("Closed webserver on http://127.0.0.1:5000/")
