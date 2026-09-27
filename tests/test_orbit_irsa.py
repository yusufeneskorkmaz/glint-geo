"""Synthetic geometry, age, fallback, and URL checks ported from the survey."""
from datetime import datetime, timedelta, timezone

import pytest

from glint_geo.irsa import angular_separation, gaia_predicate, inside, product_url
from glint_geo.orbit import Element, acquire_elements, select_element

UTC = timezone.utc


def test_select_nearest_element_and_reject_stale():
    anchor = datetime(2024, 1, 1, tzinfo=UTC)
    elements = [Element(anchor, {}, "synthetic"),
                Element(anchor+timedelta(days=5), {}, "synthetic")]
    assert select_element(elements, anchor+timedelta(days=2)) == elements[0]
    with pytest.raises(ValueError, match="within 3 days"):
        select_element(elements, anchor+timedelta(days=9))


def test_historical_request_without_credentials_refuses_before_network(tmp_path, monkeypatch):
    monkeypatch.delenv("SPACETRACK_USER", raising=False)
    monkeypatch.delenv("SPACETRACK_PASSWORD", raising=False)
    now = datetime(2026, 9, 27, tzinfo=UTC)
    with pytest.raises(ValueError, match="historical dates require Space-Track"):
        acquire_elements(26871, now-timedelta(days=8), now-timedelta(days=7),
                         tmp_path, now=now)
    assert not list(tmp_path.rglob("*"))


def test_science_url_uses_exact_quadrant():
    row = dict(filefracday="20190103089097", field=390, filtercode="zr",
               ccdid=4, qid=3)
    assert product_url(row, "sciimg") == (
        "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/"
        "2019/0103/089097/ztf_20190103089097_000390_zr_c04_o_q3_sciimg.fits")
    assert product_url(row, "mskimg").endswith("_mskimg.fits")


def test_footprint_handles_ra_wrap():
    corners = (359.8, -0.2, 0.2, -0.2, 0.2, 0.2, 359.8, 0.2)
    assert inside(0.0, 0.0, corners)
    assert not inside(1.0, 0.0, corners)
    assert angular_separation(359.9, 0, 0.1, 0) == pytest.approx(0.2, abs=1e-6)


def test_gaia_query_wraps_ra():
    predicate = gaia_predicate(359.97, -10.0, 0.1)
    assert "ra >=" in predicate and "ra <=" in predicate and " OR " in predicate
