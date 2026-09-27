"""Plot-only display behavior; the underlying extracted CSV is unchanged."""
import csv

import pytest

from glint_geo.cli import _merge_masked_spans, _plot, _plot_data


def sample_rows():
    def row(t, mag, flux, error, quality="OK"):
        return {"seconds_from_shutter_open": str(t), "duration_s": "0.2",
                "background_subtracted_DN": str(flux),
                "flux_uncertainty_DN": str(error),
                "ZTF_mag_no_color": str(mag), "mag_uncertainty": "0.1",
                "quality": quality}

    return [row(0.1, 15.0, 30, 5), row(0.3, 15.2, 30, 5),
            row(0.5, 24.0, 1, 2), row(0.7, "", "", "", "BAD_PIXEL"),
            row(0.9, "", "", "", "BRIGHT_STAR")]


def test_relative_plot_uses_high_snr_median_and_excludes_low_snr_from_limits():
    data = _plot_data(sample_rows(), absolute_magnitude=False, snr_cutoff=3)
    assert data["reference_mag"] == pytest.approx(15.1)
    assert [sample["magnitude"] for sample in data["high"]] == pytest.approx([-0.1, 0.1])
    assert len(data["low"]) == 1
    assert data["low"][0]["magnitude"] > data["ylim"][1]
    assert [reason for _, _, reason in data["masked"]] == ["BAD_PIXEL", "BRIGHT_STAR"]


def test_absolute_plot_preserves_header_magnitude_and_renders(tmp_path):
    rows = sample_rows()
    data = _plot_data(rows, absolute_magnitude=True, snr_cutoff=3)
    assert data["reference_mag"] == 0
    assert data["high"][0]["magnitude"] == 15.0
    curve = tmp_path / "curve.csv"
    with curve.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    original_curve = curve.read_bytes()
    relative = tmp_path / "relative.png"
    absolute = tmp_path / "absolute.png"
    _plot(curve, relative)
    _plot(curve, absolute, absolute_magnitude=True)
    assert relative.read_bytes().startswith(b"\x89PNG")
    assert absolute.read_bytes().startswith(b"\x89PNG")
    assert relative.read_bytes() != absolute.read_bytes()
    assert curve.read_bytes() == original_curve


def test_invalid_cutoff_rejected():
    with pytest.raises(ValueError, match="SNR cutoff"):
        _plot_data(sample_rows(), absolute_magnitude=False, snr_cutoff=0)


def test_adjacent_masked_bins_share_one_shaded_interval():
    spans = [(0.0, 0.2, "BRIGHT_STAR"), (0.2, 0.4, "BRIGHT_STAR"),
             (0.4, 0.6, "BAD_PIXEL")]
    assert _merge_masked_spans(spans) == [(0.0, 0.4, "BRIGHT_STAR"),
                                         (0.4, 0.6, "BAD_PIXEL")]
