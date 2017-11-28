# encoding: utf-8
from microSALT.reports import app
import microSALT.reports.views

if __name__ == '__main__':
    app.config['DEBUG'] = True
    app.run()

