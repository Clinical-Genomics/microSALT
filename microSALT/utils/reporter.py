"""Generates various reports by tapping into Flask and mySQL
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python
import weasyprint
from multiprocessing import Process
from microSALT.server.views import app

class Reporter():

  def __init__(self, config, log):
    self.name = ""
    self.config = config
    self.logger = log
    self.server = Process(target=app.run)
    self.server.start()

  def gen_pdf(self, name):
    self.name = name
    weasyprint.HTML('http://127.0.0.1:5000/microSALT/{}/all'.format(self.name))\
    .write_pdf('{}.pdf'.format(self.name))

  def gen_excel(self):
    pass

  def kill_flask(self):
    self.server.terminate()
    self.server.join()

