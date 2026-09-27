"""Opt-in public-network smoke check; never part of routine offline tests."""
import os

import httpx
import pytest


@pytest.mark.live
@pytest.mark.skipif(os.environ.get("GLINT_GEO_LIVE") != "1", reason="set GLINT_GEO_LIVE=1")
def test_public_irsa_service_reachable():
    response = httpx.get("https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci",
                         params={"FORMAT": "METADATA", "ct": "csv"}, timeout=30)
    assert response.status_code == 200
