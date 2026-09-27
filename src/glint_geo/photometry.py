"""Registered streak photometry with per-segment quality masks."""
from __future__ import annotations

import math

import numpy as np
from scipy.ndimage import map_coordinates, median_filter

BAD_BITS = sum(1 << bit for bit in (*range(2, 11), 12))


def field_star_calibration_status(header) -> dict:
    """Audit ZTF's header field-star solution before reporting magnitudes."""
    zp = float(header.get("MAGZP", math.nan))
    n = int(header.get("NMATCHES", 0))
    rms = float(header.get("MAGZPRMS", math.nan))
    apcor = float(header.get("APCOR6", math.nan))
    return {"supported": bool(math.isfinite(zp) and n >= 10 and
                              math.isfinite(rms) and math.isfinite(apcor) and
                              abs(apcor) <= 0.05),
            "magzp": zp, "n_field_stars": n, "zp_rms_mag": rms,
            "large_aperture_correction_mag": apcor}


def ztf_magnitude(total_dn: float, magzp: float) -> float | None:
    if total_dn <= 0 or not math.isfinite(total_dn) or not math.isfinite(magzp):
        return None
    return magzp - 2.5 * math.log10(total_dn)


def _line_samples(image: np.ndarray, start: np.ndarray, end: np.ndarray,
                  shift: np.ndarray, halfwidth: int) -> tuple[np.ndarray, np.ndarray]:
    length = float(np.linalg.norm(end - start))
    n = max(2, int(round(length)))
    direction = (end - start) / length
    normal = np.array([-direction[1], direction[0]])
    base = start[:, None] + direction[:, None] * np.arange(n)[None, :] + shift[:, None]
    aperture = np.arange(-halfwidth, halfwidth + 1)
    side = np.r_[np.arange(halfwidth + 6, halfwidth + 11),
                 -np.arange(halfwidth + 6, halfwidth + 11)]

    def take(widths):
        points = base[:, :, None] + normal[:, None, None] * widths[None, None, :]
        return map_coordinates(image, [points[1].ravel(), points[0].ravel()],
                               order=1, mode="constant", cval=np.nan).reshape(n, len(widths))

    a, b = take(aperture), take(side)
    valid = np.all(np.isfinite(a), axis=1) & np.all(np.isfinite(b), axis=1)
    flux, scatter = np.full(n, np.nan), np.full(n, np.nan)
    if np.any(valid):
        sky = np.median(b[valid], axis=1)
        flux[valid] = np.sum(a[valid] - sky[:, None], axis=1)
        scatter[valid] = np.std(b[valid], axis=1)
    return flux, scatter


def trace_profile(image: np.ndarray, start: np.ndarray, end: np.ndarray,
                  *, seeing_px: float, search_px: int = 25,
                  gain_e_per_dn: float = 1.0) -> dict:
    """Register a coherent streak near the prediction, including shutter edges."""
    image = np.asarray(image, dtype=float)
    start, end = np.asarray(start, dtype=float), np.asarray(end, dtype=float)
    if image.ndim != 2 or not np.all(np.isfinite(start)) or not np.all(np.isfinite(end)):
        raise ValueError("invalid cutout or track")
    length = float(np.linalg.norm(end-start))
    if length < 10:
        raise ValueError("predicted streak too short")
    hw = max(2, int(math.ceil(seeing_px * 0.9)))
    diff = np.diff(image, axis=0)
    finite = diff[np.isfinite(diff)]
    if not len(finite):
        raise ValueError("cutout has no finite pixels")
    sigma_px = max(1e-6, 1.4826 * np.median(np.abs(finite - np.median(finite))) / math.sqrt(2))
    direction = (end-start)/length
    normal = np.array([-direction[1], direction[0]])

    def score(shift):
        flux, _ = _line_samples(image, start, end, shift, hw)
        samples = flux[np.isfinite(flux)]
        if len(samples) < 0.8 * len(flux):
            return -math.inf
        lo, hi = np.quantile(samples, [0.1, 0.9])
        return float(np.mean(samples[(samples >= lo) & (samples <= hi)]))

    _, cross = max((score(normal*offset), offset) for offset in range(-search_px, search_px+1))
    shift = normal * cross
    pad = min(150, max(120, int(length * 0.35)))
    extended, _ = _line_samples(image, start-direction*pad, end+direction*pad, shift, hw)
    smooth = median_filter(np.nan_to_num(extended, nan=0), size=7, mode="nearest")
    threshold = max(4*sigma_px*math.sqrt(2*hw+1), 0.08*float(np.max(smooth)))
    bounds = np.flatnonzero(np.diff(np.r_[False, smooth > threshold, False].astype(int)))
    runs = [(int(a), int(b)) for a, b in bounds.reshape(-1, 2)]
    long_runs = [run for run in runs if run[1] - run[0] >= 0.55*length]
    along = max(long_runs, key=lambda r: r[1]-r[0])[0]-pad if long_runs else 0
    shift = shift + direction*along
    flux, _ = _line_samples(image, start, end, shift, hw)
    good = np.isfinite(flux)
    if not np.any(good):
        return {"present": False, "shift": shift, "flux": flux,
                "pixel_sigma": sigma_px, "aperture_halfwidth_px": hw,
                "snr": 0.0, "cross_track_offset_px": int(cross),
                "along_track_offset_px": int(along)}
    lo, hi = np.quantile(flux[good], [0.1, 0.9])
    core = flux[good & (flux >= lo) & (flux <= hi)]
    total = float(np.mean(core) * len(core))
    variance = sigma_px**2 * (2*hw+1) * len(core) * 1.7 + max(total, 0)/max(gain_e_per_dn, 1e-6)
    snr = total/math.sqrt(variance)
    return {"present": bool(snr >= 10 and np.mean(flux[good] > 0) >= 0.65),
            "shift": shift, "flux": flux, "pixel_sigma": sigma_px,
            "aperture_halfwidth_px": hw, "snr": snr,
            "cross_track_offset_px": int(cross), "along_track_offset_px": int(along)}


def _sample_strip(image: np.ndarray, mask: np.ndarray, start: np.ndarray,
                  end: np.ndarray, shift: np.ndarray, halfwidth: int,
                  stars_xy_g: np.ndarray, star_radius_px: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    length = float(np.linalg.norm(end - start))
    n = int(round(length))
    direction = (end - start) / length
    normal = np.array([-direction[1], direction[0]])
    base = start[:, None] + shift[:, None] + direction[:, None] * np.arange(n)[None, :]
    aperture = np.arange(-halfwidth, halfwidth + 1)
    side = np.r_[np.arange(halfwidth+6, halfwidth+11), -np.arange(halfwidth+6, halfwidth+11)]

    def take(widths, data, order, cval):
        points = base[:, :, None] + normal[:, None, None] * widths[None, None, :]
        return map_coordinates(data, [points[1].ravel(), points[0].ravel()],
                               order=order, mode="constant", cval=cval).reshape(n, len(widths))

    a, b = take(aperture, image, 1, np.nan), take(side, image, 1, np.nan)
    ma, mb = take(aperture, mask, 0, BAD_BITS), take(side, mask, 0, BAD_BITS)
    bad = np.any((ma.astype(np.uint16) & BAD_BITS) != 0, axis=1)
    bad |= np.any((mb.astype(np.uint16) & BAD_BITS) != 0, axis=1)
    bad |= ~np.all(np.isfinite(a), axis=1) | ~np.all(np.isfinite(b), axis=1)
    star_bad = np.zeros(n, dtype=bool)
    for x, y, mag in stars_xy_g:
        if mag >= 16:
            continue
        offset = np.array([x, y]) - (start + shift)
        along = float(offset @ direction)
        across = float(offset @ normal)
        if abs(across) <= star_radius_px:
            star_bad |= np.abs(np.arange(n)-along) <= star_radius_px
    sky = np.median(b, axis=1)
    flux = np.sum(a-sky[:, None], axis=1)
    scatter = np.std(b, axis=1)
    flux[bad | star_bad] = np.nan
    scatter[bad | star_bad] = np.nan
    reason = np.where(bad, "BAD_PIXEL", np.where(star_bad, "BRIGHT_STAR", "OK"))
    return flux, scatter, reason


def extract_curve(image: np.ndarray, mask: np.ndarray, start: np.ndarray,
                  end: np.ndarray, header, stars_xy_g: np.ndarray,
                  *, exptime_s: float, seeing_arcsec: float | None) -> tuple[list[dict], dict]:
    """Return sub-second bins and a calibration/quality summary."""
    image = np.asarray(image, dtype=float)
    mask = np.asarray(mask, dtype=np.uint16)
    if image.shape != mask.shape or image.ndim != 2:
        raise ValueError("science and mask cutout dimensions differ")
    calibration = field_star_calibration_status(header)
    if not calibration["supported"]:
        raise ValueError("image header lacks a supported field-star zero point")
    gain = float(header.get("GAIN", 1.0))
    seeing_px = (seeing_arcsec or 2.0)/1.01
    profile = trace_profile(image, start, end, seeing_px=seeing_px,
                            gain_e_per_dn=gain)
    if not profile["present"]:
        raise ValueError("no coherent streak detected at the predicted location")
    shift = profile["shift"]
    margin = min(start[0]+shift[0], start[1]+shift[1], end[0]+shift[0], end[1]+shift[1],
                 image.shape[1]-1-start[0]-shift[0], image.shape[0]-1-start[1]-shift[1],
                 image.shape[1]-1-end[0]-shift[0], image.shape[0]-1-end[1]-shift[1])
    if margin < 5:
        raise ValueError("registered streak is truncated by cutout edge")
    stars_xy_g = np.asarray(stars_xy_g, dtype=float).reshape(-1, 3)
    flux, scatter, reasons = _sample_strip(image, mask, start, end, shift,
                                            profile["aperture_halfwidth_px"],
                                            stars_xy_g,
                                            profile["aperture_halfwidth_px"] + 11 + seeing_px)
    n = len(flux)
    pixels_per_second = n/exptime_s
    bin_pixels = max(1, min(3, int(math.floor(0.9*pixels_per_second))))
    cadence = bin_pixels/pixels_per_second
    if cadence >= 1:
        raise ValueError("predicted trail cannot support sub-second bins")
    rows = []
    aperture_pixels = 2*profile["aperture_halfwidth_px"]+1
    for left in range(0, n, bin_pixels):
        right = min(n, left+bin_pixels)
        segment = flux[left:right]
        statuses = reasons[left:right]
        duration = (right-left)/pixels_per_second
        t = (left+right)/(2*pixels_per_second)
        if np.any(statuses != "OK") or not np.all(np.isfinite(segment)):
            reason = "BAD_PIXEL" if "BAD_PIXEL" in statuses else "BRIGHT_STAR"
            rows.append({"seconds_from_shutter_open": t, "duration_s": duration,
                         "background_subtracted_DN": None, "flux_uncertainty_DN": None,
                         "ZTF_mag_no_color": None, "mag_uncertainty": None,
                         "quality": reason})
            continue
        value = float(np.sum(segment))
        sigma = max(profile["pixel_sigma"], float(np.nanmedian(scatter[left:right])))
        error = math.sqrt(sigma**2 * aperture_pixels * len(segment) * 1.7 +
                          max(value, 0)/max(gain, 1e-6))
        mag = ztf_magnitude(value * exptime_s/duration, calibration["magzp"])
        rows.append({"seconds_from_shutter_open": t, "duration_s": duration,
                     "background_subtracted_DN": value,
                     "flux_uncertainty_DN": error,
                     "ZTF_mag_no_color": mag,
                     "mag_uncertainty": 1.085736*error/value if value > 0 else None,
                     "quality": "OK"})
    good = sum(row["quality"] == "OK" for row in rows)
    if good == 0:
        raise ValueError("all streak segments are masked")
    summary = {**calibration, "gain_e_per_dn": gain, "integrated_snr": profile["snr"],
               "cross_track_offset_px": profile["cross_track_offset_px"],
               "along_track_offset_px": profile["along_track_offset_px"],
               "cadence_s": cadence, "n_bins": len(rows), "n_good_bins": good,
               "n_bad_pixel_bins": sum(row["quality"] == "BAD_PIXEL" for row in rows),
               "n_bright_star_bins": sum(row["quality"] == "BRIGHT_STAR" for row in rows)}
    return rows, summary
