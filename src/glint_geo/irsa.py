"""Public IRSA metadata, footprint matching, cutouts, and star catalogs."""
from __future__ import annotations

import csv
from datetime import datetime, timedelta, timezone
from hashlib import sha256
import io
import math
from pathlib import Path

import httpx
from skyfield.api import load

from .orbit import radec, select_element, utc

IRSA = "https://irsa.ipac.caltech.edu"
JD_UNIX = 2440587.5
CORNER_KEYS = ("ra1", "dec1", "ra2", "dec2", "ra3", "dec3", "ra4", "dec4")


def jd(value: datetime) -> float:
    return JD_UNIX + utc(value).timestamp() / 86400.0


def from_jd(value: float) -> datetime:
    return datetime.fromtimestamp((value - JD_UNIX) * 86400.0, timezone.utc)


def angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    a, b, c, d = map(math.radians, (ra1, dec1, ra2, dec2))
    cos_sep = math.sin(b) * math.sin(d) + math.cos(b) * math.cos(d) * math.cos(a-c)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_sep))))


def inside(ra: float, dec: float, corners: tuple[float, ...]) -> bool:
    """Convex quadrant containment with right-ascension wrapping."""
    x = [((r - ra + 180.0) % 360.0 - 180.0) * math.cos(math.radians(dec))
         for r in corners[::2]]
    y = [d - dec for d in corners[1::2]]
    cross = [x[i] * y[(i+1) % 4] - y[i] * x[(i+1) % 4] for i in range(4)]
    return min(cross) >= -1e-10 or max(cross) <= 1e-10


def _tap(client: httpx.Client, query: str) -> list[dict[str, str]]:
    try:
        response = client.get("/TAP/sync", params={"QUERY": query, "FORMAT": "csv"})
    except httpx.HTTPError:
        raise RuntimeError("IRSA TAP connection failed") from None
    if response.status_code != 200 or response.is_redirect or len(response.content) > 64 * 1024 * 1024:
        raise RuntimeError("IRSA TAP request failed or exceeded 64 MiB")
    text = response.content.decode("utf-8-sig")
    if text.lstrip().startswith("<"):
        raise RuntimeError("IRSA TAP returned an error instead of CSV")
    return list(csv.DictReader(io.StringIO(text)))


def _raw_query(lo: datetime, hi: datetime, exposure_id: int | None) -> str:
    filter_id = f" AND expid={exposure_id}" if exposure_id is not None else ""
    return ("SELECT DISTINCT expid,obsjd,telra,teldec,exptime "
            "FROM ztf.ztf_current_meta_raw "
            f"WHERE obsjd>={jd(lo):.9f} AND obsjd<{jd(hi):.9f} "
            "AND imgtypecode='o' AND ipac_pub_date IS NOT NULL" + filter_id)


def _science_query(ids: list[int]) -> str:
    if not ids:
        raise ValueError("empty exposure list")
    return ("SELECT expid,ccdid,qid,obsjd,exptime,ra1,dec1,ra2,dec2,"
            "ra3,dec3,ra4,dec4,field,filtercode,filefracday,seeing,infobits "
            "FROM ztf.ztf_current_meta_sci WHERE expid IN (" +
            ",".join(map(str, ids)) + ")")


def find_exposures(identifier: int, start: datetime, end: datetime,
                   elements: list, *, exposure_id: int | None = None,
                   limit: int | None = None) -> list[dict]:
    """Find public exact-quadrant intersections; reject stale elements per epoch."""
    start, end = utc(start), utc(end)
    if end < start:
        raise ValueError("end precedes start")
    ts = load.timescale(builtin=True)
    coarse: dict[int, dict] = {}
    with httpx.Client(base_url=IRSA, follow_redirects=False, trust_env=False,
                      timeout=180.0) as client:
        day = start.replace(hour=0, minute=0, second=0, microsecond=0)
        while day <= end:
            for row in _tap(client, _raw_query(day, day + timedelta(days=1), exposure_id)):
                expid = int(row["expid"])
                obsjd = float(row["obsjd"])
                exptime = float(row["exptime"])
                middle = from_jd(obsjd + exptime / 172800.0)
                if not (start <= middle <= end):
                    continue
                try:
                    element = select_element(elements, middle)
                except ValueError:
                    continue
                ra, dec = radec(element, middle, ts)
                if angular_separation(ra, dec, float(row["telra"]), float(row["teldec"])) <= 6.5:
                    coarse[expid] = {"expid": expid, "obsjd": obsjd,
                                     "exptime_s": exptime, "midpoint_utc": middle,
                                     "ra_deg": ra, "dec_deg": dec,
                                     "element": element}
            day += timedelta(days=1)
        matches = []
        ids = sorted(coarse)
        for offset in range(0, len(ids), 200):
            for row in _tap(client, _science_query(ids[offset:offset+200])):
                base = coarse.get(int(row["expid"]))
                if base is None:
                    raise RuntimeError("IRSA returned an unrequested exposure ID")
                corners = tuple(float(row[k]) for k in CORNER_KEYS)
                if not inside(base["ra_deg"], base["dec_deg"], corners):
                    continue
                shutter = from_jd(base["obsjd"])
                close = shutter + timedelta(seconds=base["exptime_s"])
                try:
                    first = radec(select_element(elements, shutter), shutter, ts)
                    last = radec(select_element(elements, close), close, ts)
                except ValueError:
                    continue
                if not (inside(*first, corners) and inside(*last, corners)):
                    continue
                matches.append({
                    "norad_id": identifier, "expid": base["expid"],
                    "ccdid": int(row["ccdid"]), "qid": int(row["qid"]),
                    "field": int(row["field"]), "filtercode": row["filtercode"],
                    "filefracday": row["filefracday"],
                    "obsjd": base["obsjd"], "exptime_s": base["exptime_s"],
                    "observation_utc": shutter.isoformat(),
                    "midpoint_ra_deg": base["ra_deg"], "midpoint_dec_deg": base["dec_deg"],
                    "start_radec_deg": [*first], "end_radec_deg": [*last],
                    "element_source": base["element"].source,
                    "element_age_days": abs((base["element"].epoch - base["midpoint_utc"]).total_seconds())/86400,
                    "seeing_arcsec": float(row["seeing"]) if row.get("seeing") else None,
                    "infobits": int(row["infobits"]) if row.get("infobits") else None,
                })
    matches.sort(key=lambda r: (r["obsjd"], r["expid"], r["ccdid"], r["qid"]))
    return matches[:limit] if limit is not None else matches


def product_url(row: dict, kind: str) -> str:
    if kind not in {"sciimg", "mskimg"}:
        raise ValueError("only science and mask images are supported")
    day = str(row["filefracday"])
    if len(day) != 14 or not day.isascii() or not day.isdigit():
        raise ValueError("invalid IRSA filefracday")
    name = (f"ztf_{day}_{int(row['field']):06d}_{row['filtercode']}_"
            f"c{int(row['ccdid']):02d}_o_q{int(row['qid'])}_{kind}.fits")
    return f"{IRSA}/ibe/data/ztf/products/sci/{day[:4]}/{day[4:8]}/{day[8:]}/{name}"


def _cache_bytes(path: Path, payload: bytes) -> str:
    digest = sha256(payload).hexdigest()
    path.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
    if path.exists():
        if sha256(path.read_bytes()).hexdigest() != digest:
            raise RuntimeError("cached response differs; refusing overwrite")
    else:
        path.write_bytes(payload)
    return digest


def cutout_size_arcsec(row: dict) -> int:
    a, b = row["start_radec_deg"], row["end_radec_deg"]
    length = angular_separation(*a, *b) * 3600
    size = max(600, math.ceil(length + 200))
    if size > 1600:
        raise ValueError("predicted trail is too large for a safe cutout")
    return size


def fetch_cutouts(row: dict, cache: Path) -> dict:
    """Fetch two public FITS cutouts, recording URL and SHA-256 for each."""
    size = cutout_size_arcsec(row)
    center = f"{row['midpoint_ra_deg']:.8f},{row['midpoint_dec_deg']:.8f}"
    key = f"{row['norad_id']}_{row['expid']}_c{row['ccdid']:02d}q{row['qid']}"
    products = {}
    with httpx.Client(follow_redirects=False, trust_env=False, timeout=180.0) as client:
        for kind in ("sciimg", "mskimg"):
            rel = f"cutouts/{key}_{kind}.fits"
            path = cache / rel
            url = product_url(row, kind)
            params = {"center": center, "size": f"{size}arcsec", "gzip": "false"}
            request_url = str(httpx.URL(url, params=params))
            if path.exists():
                payload = path.read_bytes()
            else:
                try:
                    response = client.get(url, params=params)
                except httpx.HTTPError:
                    raise RuntimeError("IRSA cutout connection failed") from None
                if response.status_code != 200 or response.is_redirect:
                    raise RuntimeError(f"public IRSA {kind} cutout unavailable")
                payload = response.content
            if not payload.startswith(b"SIMPLE") or len(payload) > 128 * 1024 * 1024:
                raise RuntimeError("IRSA returned an invalid or oversized FITS cutout")
            products[kind] = {"url": request_url, "sha256": _cache_bytes(path, payload),
                              "cache": rel}
    return products


def gaia_predicate(ra_deg: float, dec_deg: float, radius_deg: float) -> str:
    half_ra = radius_deg / max(math.cos(math.radians(dec_deg)), 0.2)
    lo, hi = ra_deg - half_ra, ra_deg + half_ra
    if lo < 0:
        ra_part = f"(ra >= {lo+360:.8f} OR ra <= {hi:.8f})"
    elif hi >= 360:
        ra_part = f"(ra >= {lo:.8f} OR ra <= {hi-360:.8f})"
    else:
        ra_part = f"ra BETWEEN {lo:.8f} AND {hi:.8f}"
    return f"{ra_part} AND dec BETWEEN {dec_deg-radius_deg:.8f} AND {dec_deg+radius_deg:.8f}"


def fetch_bright_stars(row: dict, cache: Path) -> dict:
    """Cache a public IRSA-hosted Gaia bright-star catalog for segment masking."""
    key = f"{row['norad_id']}_{row['expid']}_c{row['ccdid']:02d}q{row['qid']}"
    rel = f"catalogs/{key}_gaia.csv"
    path = cache / rel
    query = ("SELECT ra,dec,phot_g_mean_mag FROM gaia_dr3_source WHERE " +
             gaia_predicate(row["midpoint_ra_deg"], row["midpoint_dec_deg"],
                            cutout_size_arcsec(row)/7200 + 0.02) +
             " AND phot_g_mean_mag < 16")
    if path.exists():
        payload = path.read_bytes()
    else:
        with httpx.Client(base_url=IRSA, follow_redirects=False,
                          trust_env=False, timeout=120.0) as client:
            records = _tap(client, query)
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fieldnames=["ra", "dec", "phot_g_mean_mag"])
        writer.writeheader()
        writer.writerows(records)
        payload = stream.getvalue().encode()
    if not payload.startswith(b"ra,dec,phot_g_mean_mag"):
        raise RuntimeError("invalid bright-star catalog response")
    return {"cache": rel, "sha256": _cache_bytes(path, payload)}
