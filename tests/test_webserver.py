#!/usr/bin/env python

import requests
import time
from unittest.mock import patch

from flask import Flask

from microsalt.utils.reporter import Reporter
from microsalt.server.views import *


def test_webserver(reporter: Reporter):
    reporter.start_web()
    reporter.kill_flask()


def test_pages(reporter: Reporter):
    reporter.start_web()
    # Valid pages with available data
    time.sleep(3)
    a = requests.get("http://127.0.0.1:5001/", allow_redirects=True)
    assert a.status_code == 200

    time.sleep(0.15)
    b = requests.get("http://127.0.0.1:5001/microsalt/", allow_redirects=True)
    assert b.status_code == 200

    time.sleep(0.15)
    c = requests.get("http://127.0.0.1:5001/microsalt/AAA1234", allow_redirects=True)
    assert c.status_code == 200

    time.sleep(0.15)
    e = requests.get("http://127.0.0.1:5001/microsalt/AAA1234/typing/all", allow_redirects=True)
    assert e.status_code in [200, 500]

    time.sleep(0.15)
    d = requests.get("http://127.0.0.1:5001/microsalt/AAA1234/qc", allow_redirects=True)
    assert d.status_code in [200, 500]

    # Valid pages with unavailable data
    time.sleep(0.15)
    f = requests.get(
        "http://127.0.0.1:5001/microsalt/AAA1234/typing/escherichia_coli", allow_redirects=True
    )
    assert f.status_code in [200, 500]

    time.sleep(0.15)
    g = requests.get("http://127.0.0.1:5001/microsalt/STtracker/all", allow_redirects=True)
    assert g.status_code in [200, 500]

    time.sleep(0.15)
    h = requests.get("http://127.0.0.1:5001/microsalt/STtracker/cust000", allow_redirects=True)
    assert h.status_code in [200, 500]


@patch("microsalt.server.views.render_template")
def test_index_views(renderpatch, app: Flask):
    renderpatch.return_value = "ok"
    with app.app_context():
        start = start_page()
        assert start == "ok"
        reroute = reroute_page()
        assert reroute == "ok"


@patch("microsalt.server.views.render_template")
def test_project_views(renderpatch, app: Flask, dbm):
    renderpatch.return_value = "ok"
    with app.app_context():
        a = project_page("AAA1234")
        assert a == "ok"
        b = alignment_page("AAA1234")
        assert b == "ok"
        c = typing_page("AAA1234", "all")
        assert c == "ok"


@patch("microsalt.server.views.gen_add_info")
@patch("microsalt.server.views.render_template")
def test_tracker_view(renderpatch, addinfo):
    renderpatch.return_value = "ok"
    a = STtracker_page("cust000")
    assert a == "ok"
