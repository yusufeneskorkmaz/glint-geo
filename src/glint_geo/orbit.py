"""Private, local GP element acquisition and SGP4 sky prediction."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import json
import os
from pathlib import Path

import httpx
from skyfield.api import EarthSatellite, load, wgs84

UTC = timezone.utc
MAX_AGE = timedelta(days=3)
PALOMAR = wgs84.latlon(33.357336, -116.859780, elevation_m=1702.0)


def utc(value: str | datetime) -> datetime:
    result = datetime.fromisoformat(value.replace("Z", "+00:00")) if isinstance(value, str) else value
    if result.tzinfo is None:
        raise ValueError("timestamps must include a UTC offset")
    return result.astimezone(UTC)


def _provider_epoch(value: str) -> datetime:
    """Both GP services may omit a timezone suffix but specify UTC epochs."""
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return parsed.replace(tzinfo=UTC) if parsed.tzinfo is None else parsed.astimezone(UTC)


@dataclass(frozen=True)
class Element:
    epoch: datetime
    record: dict
    source: str


def select_element(elements: list[Element], when: datetime) -> Element:
    if not elements:
        raise ValueError("no orbital elements were returned for this NORAD ID")
    when = utc(when)
    best = min(elements, key=lambda e: (abs(e.epoch - when), e.epoch))
    if abs(best.epoch - when) > MAX_AGE:
        raise ValueError("no orbital element is within 3 days of the observation epoch")
    return best


def _save_private(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    if path.parent.is_symlink():
        raise ValueError("unsafe element cache directory")
    data = (json.dumps(records, sort_keys=True, separators=(",", ":")) + "\n").encode()
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(fd, "wb") as stream:
        stream.write(data)


def _space_track(identifier: int, start: datetime, end: datetime) -> list[dict]:
    username = os.environ.get("SPACETRACK_USER")
    password = os.environ.get("SPACETRACK_PASSWORD")
    if not username or not password:
        raise ValueError("both Space-Track environment variables are required")
    lo = (start - MAX_AGE).strftime("%Y-%m-%dT%H:%M:%S")
    hi = (end + MAX_AGE).strftime("%Y-%m-%dT%H:%M:%S")
    path = (f"/basicspacedata/query/class/gp_history/NORAD_CAT_ID/{identifier}/"
            f"EPOCH/{lo}--{hi}/orderby/epoch/format/json")
    try:
        with httpx.Client(base_url="https://www.space-track.org", follow_redirects=False,
                          trust_env=False, timeout=120.0) as client:
            login = client.post("/ajaxauth/login", data={"identity": username, "password": password})
            if login.status_code != 200 or login.is_redirect:
                raise RuntimeError("Space-Track authentication failed")
            response = client.get(path)
            if response.status_code != 200 or response.is_redirect:
                raise RuntimeError("Space-Track GP_HISTORY request failed")
            if len(response.content) > 64 * 1024 * 1024:
                raise RuntimeError("Space-Track response exceeds 64 MiB")
            records = response.json()
    except httpx.HTTPError as exc:
        raise RuntimeError("Space-Track connection failed") from None
    if not isinstance(records, list):
        raise RuntimeError("Space-Track returned an invalid GP_HISTORY response")
    return records


def _celestrak(identifier: int, start: datetime, end: datetime, now: datetime) -> list[dict]:
    if start < now - MAX_AGE or end > now:
        raise ValueError("historical dates require Space-Track GP_HISTORY access; CelesTrak supplies current elements only")
    try:
        with httpx.Client(follow_redirects=False, trust_env=False, timeout=45.0) as client:
            response = client.get("https://celestrak.org/NORAD/elements/gp.php",
                                  params={"CATNR": identifier, "FORMAT": "JSON"})
            if response.status_code != 200 or response.is_redirect:
                raise RuntimeError("CelesTrak current-GP request failed")
            if len(response.content) > 1024 * 1024:
                raise RuntimeError("CelesTrak response exceeds 1 MiB")
            records = response.json()
    except httpx.HTTPError:
        raise RuntimeError("CelesTrak connection failed") from None
    if not isinstance(records, list):
        raise RuntimeError("CelesTrak returned an invalid GP response")
    return records


def acquire_elements(identifier: int, start: datetime, end: datetime,
                     cache: Path, *, now: datetime | None = None) -> list[Element]:
    """Cache only inside an ignored private folder; never expose records in output manifests."""
    if not 1 <= identifier <= 999_999_999:
        raise ValueError("NORAD ID must be a positive catalog number")
    start, end = utc(start), utc(end)
    if end < start:
        raise ValueError("end precedes start")
    now = utc(now or datetime.now(UTC))
    historical = bool(os.environ.get("SPACETRACK_USER") and os.environ.get("SPACETRACK_PASSWORD"))
    source = "Space-Track GP_HISTORY" if historical else "CelesTrak current GP"
    key = f"{identifier}_{start:%Y%m%dT%H%M%S}_{end:%Y%m%dT%H%M%S}_{'history' if historical else 'current'}.json"
    path = cache / "elements" / key
    if path.exists():
        records = json.loads(path.read_text())
    else:
        records = (_space_track(identifier, start, end) if historical else
                   _celestrak(identifier, start, end, now))
        _save_private(path, records)
    elements = []
    for record in records:
        if not isinstance(record, dict) or int(record.get("NORAD_CAT_ID", -1)) != identifier:
            raise RuntimeError("element response contains an unexpected catalog ID")
        epoch = _provider_epoch(record["EPOCH"])
        elements.append(Element(epoch, record, source))
    return sorted(elements, key=lambda e: e.epoch)


def satellite(element: Element, timescale) -> EarthSatellite:
    record = element.record
    if record.get("TLE_LINE1") and record.get("TLE_LINE2"):
        return EarthSatellite(record["TLE_LINE1"], record["TLE_LINE2"],
                              str(record.get("OBJECT_NAME", "")), timescale)
    return EarthSatellite.from_omm(timescale, record)


def radec(element: Element, when: datetime, timescale=None) -> tuple[float, float]:
    timescale = timescale or load.timescale(builtin=True)
    target = satellite(element, timescale)
    ra, dec, _ = (target - PALOMAR).at(timescale.from_datetime(utc(when))).radec()
    return float(ra.degrees), float(dec.degrees)
