"""Synthetic streak tests adapted from the GEO survey image checks."""
import math

import numpy as np
import pytest

from glint_geo.photometry import (
    _line_samples, extract_curve, field_star_calibration_status,
    trace_profile, ztf_magnitude,
)


def synthetic_streak():
    image = np.full((120, 140), 100.0)
    image += np.random.default_rng(3).normal(0, 2, image.shape)
    for x in range(20, 101):
        for y in range(40, 47):
            image[y, x] += 30*np.exp(-0.5*((y-43)/1.2)**2)
    return image


def test_registers_offset_from_synthetic_streak():
    found = trace_profile(synthetic_streak(), np.array([20., 40.]),
                          np.array([100., 40.]), seeing_px=2.8, search_px=8)
    assert found["present"]
    assert abs(found["shift"][1]-3) <= 1
    assert found["snr"] > 20


def test_blank_and_bright_point_are_not_streaks():
    blank = 100+np.random.default_rng(4).normal(0, 2, (120, 120))
    point_y, point_x = np.mgrid[:120, :120]
    point = 100+10000*np.exp(-((point_x-60)**2+(point_y-40)**2)/(2*1.5**2))
    point += np.random.default_rng(17).normal(0, 2, point.shape)
    for image in (blank, point):
        assert not trace_profile(image, np.array([20., 40.]),
                                 np.array([100., 40.]), seeing_px=2.8,
                                 search_px=8)["present"]


def test_outside_cutout_is_missing():
    flux, scatter = _line_samples(np.ones((50, 50)), np.array([60., 20.]),
                                  np.array([80., 20.]), np.zeros(2), 2)
    assert np.all(np.isnan(flux))
    assert np.all(np.isnan(scatter))


def test_header_field_star_calibration_and_zero_point():
    header = {"MAGZP": 25., "NMATCHES": 223, "MAGZPRMS": .035, "APCOR6": -.004}
    assert field_star_calibration_status(header)["supported"]
    assert ztf_magnitude(10000, 25) == pytest.approx(15)
    header["APCOR6"] = -0.4
    assert not field_star_calibration_status(header)["supported"]


def test_segment_masking_preserves_other_bins():
    image = synthetic_streak()
    mask = np.zeros(image.shape, dtype=np.uint16)
    mask[43, 40:44] |= 1 << 8
    header = {"MAGZP": 25., "NMATCHES": 100, "MAGZPRMS": .04,
              "APCOR6": .0, "GAIN": 6.}
    stars = np.array([[70., 43., 14.]])
    rows, summary = extract_curve(image, mask, np.array([20., 40.]),
                                  np.array([100., 40.]), header, stars,
                                  exptime_s=5.0, seeing_arcsec=2.8)
    assert summary["cadence_s"] < 1
    assert summary["n_good_bins"] > 0
    assert summary["n_bad_pixel_bins"] > 0
    assert summary["n_bright_star_bins"] > 0
    assert any(r["quality"] == "OK" and math.isfinite(r["ZTF_mag_no_color"])
               for r in rows)
    assert all(r["background_subtracted_DN"] is None for r in rows
               if r["quality"] != "OK")


def test_unsupported_header_refuses_calibrated_curve():
    with pytest.raises(ValueError, match="field-star zero point"):
        extract_curve(synthetic_streak(), np.zeros((120, 140), dtype=np.uint16),
                      np.array([20., 40.]), np.array([100., 40.]), {},
                      np.empty((0, 3)), exptime_s=5, seeing_arcsec=2.8)
