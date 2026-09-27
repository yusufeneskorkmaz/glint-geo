"""Offline FITS/WCS integration using synthetic pixels and catalog rows."""
import csv
from hashlib import sha256
import json

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from glint_geo.pipeline import extract


def test_extract_writes_curve_and_provenance_without_orbital_records(tmp_path):
    cache = tmp_path / "cache"
    out = tmp_path / "out"
    (cache / "cutouts").mkdir(parents=True)
    (cache / "catalogs").mkdir(parents=True)
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [60, 60]
    wcs.wcs.crval = [15, 0]
    wcs.wcs.cdelt = [-1/3600, 1/3600]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    header = wcs.to_header()
    header["MAGZP"] = 25.0
    header["NMATCHES"] = 100
    header["MAGZPRMS"] = 0.04
    header["APCOR6"] = 0.0
    header["GAIN"] = 6.0
    science = 100 + np.random.default_rng(3).normal(0, 2, (120, 140))
    for x in range(20, 101):
        for y in range(37, 44):
            science[y, x] += 30*np.exp(-0.5*((y-40)/1.2)**2)
    mask = np.zeros(science.shape, dtype=np.uint16)
    mask[40, 40:44] |= 1 << 8
    sci_path = cache / "cutouts" / "science.fits"
    msk_path = cache / "cutouts" / "mask.fits"
    fits.writeto(sci_path, science, header)
    fits.writeto(msk_path, mask, header)
    catalog_path = cache / "catalogs" / "stars.csv"
    catalog_path.write_text("ra,dec,phot_g_mean_mag\n")
    ra, dec = wcs.pixel_to_world_values([20, 100], [40, 40])
    item = {
        "norad_id": 26871, "expid": 1, "ccdid": 1, "qid": 1,
        "exptime_s": 5.0, "seeing_arcsec": 2.8,
        "start_radec_deg": [float(ra[0]), float(dec[0])],
        "end_radec_deg": [float(ra[1]), float(dec[1])],
        "products": {
            "sciimg": {"cache": "cutouts/science.fits", "sha256": sha256(sci_path.read_bytes()).hexdigest(), "url": "https://irsa.ipac.caltech.edu/synthetic"},
            "mskimg": {"cache": "cutouts/mask.fits", "sha256": sha256(msk_path.read_bytes()).hexdigest(), "url": "https://irsa.ipac.caltech.edu/synthetic"},
        },
        "bright_star_catalog": {"cache": "catalogs/stars.csv", "sha256": sha256(catalog_path.read_bytes()).hexdigest()},
    }
    input_manifest = tmp_path / "manifest.json"
    input_manifest.write_text(json.dumps({"format_version": 1, "streaks": [item]}))
    output_manifest = extract(input_manifest, cache, out)
    record = json.loads(output_manifest.read_text())["streaks"][0]
    assert record["status"] == "EXTRACTED"
    assert record["calibration"]["n_good_bins"] > 0
    assert record["calibration"]["n_bad_pixel_bins"] > 0
    curve = out / record["light_curve"]["file"]
    assert sha256(curve.read_bytes()).hexdigest() == record["light_curve"]["sha256"]
    rows = list(csv.DictReader(curve.open()))
    assert all(float(row["duration_s"]) < 1 for row in rows)
    assert "TLE_LINE" not in output_manifest.read_text()
