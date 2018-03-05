"""Generates various reports by tapping into Flask and mySQL
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import time
import weasyprint
from multiprocessing import Process
from microSALT.server.views import app, session
from microSALT.store.orm_models import Samples

class Reporter():

  def __init__(self, config, log):
    self.name = ""
    self.config = config
    self.logger = log
    self.server = Process(target=app.run)

  def gen_pdf(self, name):
    self.server.start()
    #Hinders requests before server goes up
    time.sleep(0.05)
    self.name = name
    weasyprint.HTML('http://127.0.0.1:5000/microSALT/{}/all'.format(self.name))\
    .write_pdf('{}.pdf'.format(self.name))
    self.kill_flask()
    self.logger.info("Created {}.pdf in current directory".format(name))

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

