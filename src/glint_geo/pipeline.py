"""Fetch, extract, and provenance output for public streaks."""
from __future__ import annotations

import csv
from hashlib import sha256
import io
import json
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from . import __version__
from .irsa import fetch_bright_stars, fetch_cutouts, find_exposures
from .orbit import acquire_elements, utc
from .photometry import extract_curve


def _json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def _safe_path(root: Path, relative: str) -> Path:
    path = (root / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValueError("manifest contains an unsafe cache path")
    return path


def fetch(identifier: int, start: str, end: str, cache: Path,
          output: Path, *, exposure_id: int | None = None,
          limit: int | None = None) -> Path:
    if limit is not None and limit < 1:
        raise ValueError("limit must be positive")
    begin, finish = utc(start), utc(end)
    elements = acquire_elements(identifier, begin, finish, cache)
    candidates = find_exposures(identifier, begin, finish, elements,
                                exposure_id=exposure_id, limit=limit)
    records = []
    for candidate in candidates:
        products = fetch_cutouts(candidate, cache)
        stars = fetch_bright_stars(candidate, cache)
        records.append({**candidate, "products": products, "bright_star_catalog": stars,
                        "status": "FETCHED"})
    path = output / "manifest.json"
    _json(path, {"format_version": 1, "software_version": __version__,
                 "norad_id": identifier, "start_utc": begin.isoformat(),
                 "end_utc": finish.isoformat(), "streaks": records})
    return path


def _verify(path: Path, expected: str) -> None:
    if not path.is_file() or sha256(path.read_bytes()).hexdigest() != expected:
        raise ValueError("cached input missing or differs from provenance hash")


def _csv_bytes(rows: list[dict]) -> bytes:
    fields = ["seconds_from_shutter_open", "duration_s", "background_subtracted_DN",
              "flux_uncertainty_DN", "ZTF_mag_no_color", "mag_uncertainty", "quality"]
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: ("" if row[key] is None else
                               f"{row[key]:.6f}" if isinstance(row[key], float) else row[key])
                         for key in fields})
    return stream.getvalue().encode()


def extract(manifest_path: Path, cache: Path, output: Path) -> Path:
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("format_version") != 1:
        raise ValueError("unsupported manifest format")
    output.mkdir(parents=True, exist_ok=True)
    for item in manifest["streaks"]:
        sci = item["products"]["sciimg"]
        msk = item["products"]["mskimg"]
        catalog = item["bright_star_catalog"]
        science_path, mask_path = _safe_path(cache, sci["cache"]), _safe_path(cache, msk["cache"])
        catalog_path = _safe_path(cache, catalog["cache"])
        for path, record in ((science_path, sci), (mask_path, msk), (catalog_path, catalog)):
            _verify(path, record["sha256"])
        with fits.open(science_path, memmap=False) as hdul:
            image = np.asarray(hdul[0].data, dtype=float)
            header = hdul[0].header.copy()
        with fits.open(mask_path, memmap=False) as hdul:
            mask = np.asarray(hdul[0].data, dtype=np.uint16)
        wcs = WCS(header).celestial
        ra = [item["start_radec_deg"][0], item["end_radec_deg"][0]]
        dec = [item["start_radec_deg"][1], item["end_radec_deg"][1]]
        x, y = wcs.world_to_pixel_values(ra, dec)
        start, end = np.array([x[0], y[0]]), np.array([x[1], y[1]])
        stars = list(csv.DictReader(io.StringIO(catalog_path.read_text())))
        if stars:
            sx, sy = wcs.world_to_pixel_values([float(s["ra"]) for s in stars],
                                                [float(s["dec"]) for s in stars])
            stars_xy_g = np.column_stack([sx, sy, [float(s["phot_g_mean_mag"]) for s in stars]])
        else:
            stars_xy_g = np.empty((0, 3))
        rows, calibration = extract_curve(image, mask, start, end, header, stars_xy_g,
                                          exptime_s=item["exptime_s"],
                                          seeing_arcsec=item["seeing_arcsec"])
        filename = (f"{item['norad_id']}_{item['expid']}_"
                    f"c{item['ccdid']:02d}q{item['qid']}.csv")
        curve_path = output / filename
        payload = _csv_bytes(rows)
        curve_path.write_bytes(payload)
        item["light_curve"] = {"file": filename, "sha256": sha256(payload).hexdigest()}
        item["calibration"] = calibration
        item["status"] = "EXTRACTED"
    result = output / "manifest.json"
    _json(result, manifest)
    return result
