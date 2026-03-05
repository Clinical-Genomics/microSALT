"""Tests for Jinja2 rendering and ORM model creation.

These tests verify:
1. ORM model instantiation with correct attributes
2. render_template produces valid HTML output
3. Each page rendering function (project_page, alignment_page, typing_page,
   STtracker_page) generates HTML that contains expected content
"""

import pytest

from microSALT.store.orm_models import (
    Collections,
    Expacs,
    Projects,
    Reports,
    Resistances,
    Samples,
    Seq_types,
    Versions,
)
from microSALT.server.views import (
    TEMPLATE_FOLDER,
    STtracker_page,
    alignment_page,
    project_page,
    render_template,
    typing_page,
)


# ---------------------------------------------------------------------------
# ORM model creation tests
# ---------------------------------------------------------------------------


class TestOrmModels:
    """Verify ORM models can be instantiated and have correct attributes."""

    def test_projects_model_attributes(self):
        p = Projects(
            CG_ID_project="TST0001",
            Customer_ID_project="CUST001",
            Customer_ID="cust001",
        )
        assert p.CG_ID_project == "TST0001"
        assert p.Customer_ID_project == "CUST001"
        assert p.Customer_ID == "cust001"

    def test_samples_model_attributes(self):
        s = Samples(
            CG_ID_sample="TST0001A1",
            CG_ID_project="TST0001",
            organism="staphylococcus_aureus",
            ST=8,
            pubmlst_ST=-1,
        )
        assert s.CG_ID_sample == "TST0001A1"
        assert s.CG_ID_project == "TST0001"
        assert s.organism == "staphylococcus_aureus"
        assert s.ST == 8

    def test_seq_types_model_attributes(self):
        st = Seq_types(
            CG_ID_sample="TST0001A1",
            loci="arcC",
            allele=6,
            contig_name="NODE_1",
            identity=100.0,
            span=1.0,
            st_predictor=True,
        )
        assert st.loci == "arcC"
        assert st.allele == 6
        assert st.identity == 100.0
        assert st.st_predictor is True

    def test_resistances_model_attributes(self):
        r = Resistances(
            CG_ID_sample="TST0001A1",
            gene="blaZ",
            instance="blaZ_1",
            contig_name="NODE_2",
            identity=99.5,
            span=0.98,
            resistance="beta-lactam",
        )
        assert r.gene == "blaZ"
        assert r.resistance == "beta-lactam"
        assert r.identity == 99.5

    def test_expacs_model_attributes(self):
        e = Expacs(
            CG_ID_sample="TST0001A1",
            gene="kpsMII",
            instance="kpsMII_1",
            contig_name="NODE_3",
            identity=98.0,
            span=0.95,
            virulence="capsule",
        )
        assert e.gene == "kpsMII"
        assert e.virulence == "capsule"

    def test_reports_model_attributes(self):
        r = Reports(
            CG_ID_project="TST0001",
            steps_aggregate="abc",
            version=1,
        )
        assert r.CG_ID_project == "TST0001"
        assert r.version == 1

    def test_versions_model_attributes(self):
        v = Versions(name="software_microSALT", version="4.3.0")
        assert v.name == "software_microSALT"
        assert v.version == "4.3.0"

    def test_collections_model_attributes(self):
        c = Collections(ID_collection="COL001", CG_ID_sample="TST0001A1")
        assert c.ID_collection == "COL001"
        assert c.CG_ID_sample == "TST0001A1"

    def test_samples_model_can_be_persisted(self, db_session):
        """Verify a model can be added to the session and queried."""
        project = Projects(
            CG_ID_project="MODEL_TEST",
            Customer_ID_project="MT001",
            Customer_ID="cust_mt",
        )
        sample = Samples(
            CG_ID_sample="MODEL_TESTA1",
            CG_ID_project="MODEL_TEST",
            organism="test_organism",
            ST=42,
            pubmlst_ST=-1,
        )
        db_session.add(project)
        db_session.add(sample)
        db_session.commit()

        result = db_session.query(Samples).filter_by(CG_ID_sample="MODEL_TESTA1").one()
        assert result.organism == "test_organism"
        assert result.ST == 42

        # Cleanup
        db_session.delete(sample)
        db_session.delete(project)
        db_session.commit()


# ---------------------------------------------------------------------------
# render_template tests
# ---------------------------------------------------------------------------


class TestRenderTemplate:
    """Verify the raw render_template function produces correct HTML."""

    def test_renders_html_string(self):
        """render_template returns a non-empty string."""
        html = render_template(
            template_folder=TEMPLATE_FOLDER,
            template_name="project_page.html",
            organisms=["all", "staphylococcus_aureus"],
            project="TST0001",
        )
        assert isinstance(html, str)
        assert len(html) > 0

    def test_renders_project_in_output(self):
        """Project ID appears in rendered project_page output."""
        html = render_template(
            template_folder=TEMPLATE_FOLDER,
            template_name="project_page.html",
            organisms=["all"],
            project="MYPROJECT",
        )
        assert "MYPROJECT" in html

    def test_renders_organisms_in_output(self):
        """Organism names appear in rendered project_page output."""
        html = render_template(
            template_folder=TEMPLATE_FOLDER,
            template_name="project_page.html",
            organisms=["all", "staphylococcus_aureus"],
            project="TST0001",
        )
        assert "staphylococcus_aureus" in html or "Staphylococcus aureus" in html

    def test_url_for_stub_injected_automatically(self):
        """url_for stub is injected when not provided, preventing TemplateError."""
        # If url_for were missing, Jinja2 would raise UndefinedError
        html = render_template(
            template_folder=TEMPLATE_FOLDER,
            template_name="project_page.html",
            organisms=["all"],
            project="TST0001",
        )
        # Should still render without error and contain HTML
        assert "<html" in html or "<!doctype" in html.lower()

    def test_render_template_uses_layout(self):
        """project_page extends layout, so layout elements appear in output."""
        html = render_template(
            template_folder=TEMPLATE_FOLDER,
            template_name="project_page.html",
            organisms=["all"],
            project="TST0001",
        )
        assert "microSALT" in html


# ---------------------------------------------------------------------------
# Page function rendering tests
# ---------------------------------------------------------------------------


class TestPageRendering:
    """Verify each view function renders real HTML with populated database."""

    def test_project_page_renders_html(self, populated_db):
        """project_page returns non-empty HTML string."""
        html = project_page("AAA1234")
        assert isinstance(html, str)
        assert len(html) > 100

    def test_project_page_contains_project_id(self, populated_db):
        """project_page HTML contains the project identifier."""
        html = project_page("AAA1234")
        assert "AAA1234" in html

    def test_project_page_contains_organisms(self, populated_db):
        """project_page lists organisms found for the project."""
        html = project_page("AAA1234")
        # Both organisms should appear (in some form)
        assert "staphylococcus_aureus" in html or "Staphylococcus aureus" in html
        assert "escherichia_coli" in html or "Escherichia coli" in html

    def test_project_page_contains_all_organism(self, populated_db):
        """project_page always includes the 'all' organism group."""
        html = project_page("AAA1234")
        assert "all" in html.lower()

    def test_alignment_page_renders_html(self, populated_db, config):
        """alignment_page returns non-empty HTML string."""
        html = alignment_page("AAA1234", threshold=config.threshold)
        assert isinstance(html, str)
        assert len(html) > 100

    def test_alignment_page_is_valid_html(self, populated_db, config):
        """alignment_page output contains basic HTML structure."""
        html = alignment_page("AAA1234", threshold=config.threshold)
        assert "<!doctype html>" in html.lower() or "<html" in html.lower()
        assert "</html>" in html.lower()

    def test_alignment_page_contains_sample_ids(self, populated_db, config):
        """alignment_page HTML contains sample identifiers."""
        html = alignment_page("AAA1234", threshold=config.threshold)
        assert "AAA1234A1" in html or "AAA1234A2" in html

    def test_typing_page_renders_html(self, populated_db, config):
        """typing_page returns non-empty HTML string."""
        html = typing_page("AAA1234", "all", threshold=config.threshold, verified_organisms=config.regex.verified_organisms)
        assert isinstance(html, str)
        assert len(html) > 100

    def test_typing_page_is_valid_html(self, populated_db, config):
        """typing_page output contains basic HTML structure."""
        html = typing_page("AAA1234", "all", threshold=config.threshold, verified_organisms=config.regex.verified_organisms)
        assert "<!doctype html>" in html.lower() or "<html" in html.lower()
        assert "</html>" in html.lower()

    def test_typing_page_contains_sample_ids(self, populated_db, config):
        """typing_page HTML contains sample identifiers."""
        html = typing_page("AAA1234", "all", threshold=config.threshold, verified_organisms=config.regex.verified_organisms)
        assert "AAA1234A1" in html or "AAA1234A2" in html

    def test_typing_page_organism_filter(self, populated_db, config):
        """typing_page filters by organism correctly."""
        html_staph = typing_page("AAA1234", "staphylococcus_aureus", threshold=config.threshold, verified_organisms=config.regex.verified_organisms)
        assert "AAA1234A1" in html_staph

    def test_sttracker_page_renders_html(self, populated_db, config):
        """STtracker_page returns non-empty HTML string."""
        html = STtracker_page("all", threshold=config.threshold)
        assert isinstance(html, str)
        assert len(html) > 100

    def test_sttracker_page_is_valid_html(self, populated_db, config):
        """STtracker_page output contains basic HTML structure."""
        html = STtracker_page("all", threshold=config.threshold)
        assert "<!doctype html>" in html.lower() or "<html" in html.lower()
        assert "</html>" in html.lower()

    def test_sttracker_page_customer_filter(self, populated_db, config):
        """STtracker_page customer='all' does not raise errors."""
        html = STtracker_page("cust000", threshold=config.threshold)
        assert isinstance(html, str)

    def test_project_page_empty_project(self, populated_db):
        """project_page for unknown project renders without error."""
        html = project_page("UNKNOWN99")
        assert isinstance(html, str)
        assert "UNKNOWN99" in html

    def test_render_alignment_page_delegates(self, populated_db, config):
        """render_alignment_page produces same output as alignment_page."""
        from microSALT.server.views import render_alignment_page

        html_a = alignment_page("AAA1234", threshold=config.threshold)
        html_b = render_alignment_page("AAA1234", threshold=config.threshold)
        assert html_a == html_b

    def test_render_typing_page_delegates(self, populated_db, config):
        """render_typing_page produces same output as typing_page."""
        from microSALT.server.views import render_typing_page

        html_a = typing_page("AAA1234", "all", threshold=config.threshold, verified_organisms=config.regex.verified_organisms)
        html_b = render_typing_page("AAA1234", "all", threshold=config.threshold, verified_organisms=config.regex.verified_organisms)
        assert html_a == html_b
