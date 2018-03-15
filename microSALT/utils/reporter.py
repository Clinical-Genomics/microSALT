"""Generates various reports by tapping into Flask and mySQL
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import requests
import time
import os
import smtplib
from os.path import basename
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.image import MIMEImage
from email.mime.application import MIMEApplication
from multiprocessing import Process
from microSALT.server.views import app, session
from microSALT.store.orm_models import Samples
from microSALT.store.lims_fetcher import LIMS_Fetcher

class Reporter():

  def __init__(self, config, log):
    self.name = ""
    self.config = config
    self.logger = log
    self.server = Process(target=app.run)

  def gen_pdf(self, name):
    self.name = name
    self.start_web()
    self.name = name
    ticketFinder = LIMS_Fetcher(self.config, self.logger)
    ticketFinder.load_lims_project_info(name)
    #wget.download("http://127.0.0.1:5000/microSALT/{}/all -o {}_microSALT.html".format(name, ticketFinder.data['Customer_ID_project']))
    r = requests.get("http://127.0.0.1:5000/microSALT/{}/all".format(name), allow_redirects=True)
    outname = "{}_microSALT.html".format(ticketFinder.data['Customer_ID_project'])
    open(outname, 'wb').write(r.content)
    self.logger.info("Generated project {} in file {}".format(name, outname))
    self.kill_flask()
    self.mail_pdf(outname)

  def mail_pdf(self, file_name):
    msg = MIMEMultipart()
    msg['Subject'] = '{} Report'.format(file_name.split('_')[0])
    msg['From'] = 'microSALT'
    msg['To'] = self.config['regex']['mail_recipient']
    #msg.attach(MIMEText(open(file_name).read()))
    part = MIMEApplication(open(file_name).read())
    part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(file_name))
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

  def gen_csv(self, name):
    excel = open("{}.csv".format(name), "w+")
    sample_info = session.query(Samples).filter(Samples.CG_ID_project==name)
    sample_info = sorted(sample_info, key=lambda sample: \
                  int(sample.CG_ID_sample.replace(sample.CG_ID_project, '')[1:]))
    excel.write("Customer_ID_sample,CG_ID_sample,organism,ST,Thresholds\n")
    for s in sample_info:
      if s.ST < 0:
        if s.ST == -1:
          s.ST = Control
        elif s.ST == -4:
          s.ST = Novel
        else:
          s.ST='None'

      threshold = 'Passed'
      for seq_type in s.seq_types:
        if seq_type.identity < 100.0:
          threshold = 'Failed'
        
      excel.write("{},{},{},{},{}\n".format(s.Customer_ID_sample, s.CG_ID_sample,\
                  s.organism.replace('_', ' ').capitalize(), s.ST,threshold))
    excel.close()
    self.logger.info("Created {}.csv in current directory".format(name))

  def kill_flask(self):
    self.server.terminate()
    self.server.join()
    self.logger.info("Closed webserver on http://127.0.0.1:5000/")

