"""Command-line interface for acquisition, extraction, and plotting."""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
from statistics import median
import sys

from .pipeline import extract, fetch


def _number(value: str | None) -> float | None:
    if not value:
        return None
    try:
        number = float(value)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def _plot_data(rows: list[dict[str, str]], *, absolute_magnitude: bool,
               snr_cutoff: float) -> dict:
    """Classify display samples without changing the extracted curve."""
    if not math.isfinite(snr_cutoff) or snr_cutoff <= 0:
        raise ValueError("SNR cutoff must be positive and finite")
    high, low, masked = [], [], []
    for row in rows:
        t, duration = _number(row.get("seconds_from_shutter_open")), _number(row.get("duration_s"))
        if t is None or duration is None or duration <= 0:
            raise ValueError("curve contains an invalid time or duration")
        quality = row.get("quality", "")
        if quality != "OK":
            masked.append((t-duration/2, t+duration/2, quality))
            continue
        flux = _number(row.get("background_subtracted_DN"))
        uncertainty = _number(row.get("flux_uncertainty_DN"))
        magnitude = _number(row.get("ZTF_mag_no_color"))
        mag_error = _number(row.get("mag_uncertainty"))
        sample = {"time": t, "magnitude": magnitude,
                  "mag_error": mag_error, "duration": duration}
        snr = flux / uncertainty if flux is not None and uncertainty and uncertainty > 0 else -math.inf
        if snr >= snr_cutoff and magnitude is not None and mag_error is not None:
            high.append(sample)
        else:
            low.append(sample)
    if not high:
        raise ValueError(f"curve has no unmasked samples with SNR >= {snr_cutoff:g}")
    reference = 0.0 if absolute_magnitude else median(s["magnitude"] for s in high)
    for sample in high + low:
        if sample["magnitude"] is not None:
            sample["magnitude"] -= reference
    limits = [sample["magnitude"] + sign * min(sample["mag_error"], 0.4)
              for sample in high for sign in (-1, 1)]
    bottom, top = min(limits), max(limits)
    padding = max(0.12, (top-bottom)*0.08)
    return {"high": high, "low": low, "masked": masked,
            "ylim": (bottom-padding, top+padding), "reference_mag": reference}


def _merge_masked_spans(spans: list[tuple[float, float, str]]) -> list[tuple[float, float, str]]:
    merged: list[tuple[float, float, str]] = []
    for left, right, reason in sorted(spans):
        if merged and merged[-1][2] == reason and left <= merged[-1][1] + 1e-6:
            previous = merged[-1]
            merged[-1] = (previous[0], max(previous[1], right), reason)
        else:
            merged.append((left, right, reason))
    return merged


def _plot(curve: Path, output: Path, *, absolute_magnitude: bool = False,
          snr_cutoff: float = 3.0) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with curve.open(newline="") as stream:
        data = _plot_data(list(csv.DictReader(stream)),
                          absolute_magnitude=absolute_magnitude,
                          snr_cutoff=snr_cutoff)
    fig, ax = plt.subplots(figsize=(9, 4))
    colors = {"BAD_PIXEL": ("#e79526", "Bad pixels"),
              "BRIGHT_STAR": ("#ad79c3", "Bright star")}
    labeled = set()
    for left, right, reason in _merge_masked_spans(data["masked"]):
        color, label = colors.get(reason, ("#999999", "Masked segment"))
        ax.axvspan(left, right, color=color, alpha=0.22,
                   label=label if label not in labeled else None, zorder=0)
        labeled.add(label)
    high = data["high"]
    ax.errorbar([s["time"] for s in high], [s["magnitude"] for s in high],
                yerr=[s["mag_error"] for s in high], fmt=".", markersize=4,
                color="#2359a0", ecolor="#9bb3cf", capsize=0,
                label=f"SNR ≥ {snr_cutoff:g}", zorder=3)
    lower, upper = data["ylim"]
    low = data["low"]
    if low:
        values = [s["magnitude"] for s in low]
        display = [min(max(value, lower), upper) if value is not None else upper
                   for value in values]
        ax.scatter([s["time"] for s in low], display, facecolors="none",
                   edgecolors="#a04b4b", marker="o", s=25, linewidths=1,
                   label=f"SNR < {snr_cutoff:g} (edge if off scale)", zorder=4)
    if not absolute_magnitude:
        ax.axhline(0, color="0.5", lw=0.7, zorder=1)
    ax.set_ylim(upper, lower)
    ax.set(xlabel="Seconds from shutter open",
           ylabel=("ZTF mag" if absolute_magnitude else "Relative ZTF mag"))
    ax.grid(alpha=0.2)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="glint-geo")
    sub = parser.add_subparsers(dest="command", required=True)
    get = sub.add_parser("fetch", help="find public ZTF intersections and cache cutouts")
    get.add_argument("--norad-id", type=int, required=True)
    get.add_argument("--start", required=True, help="timezone-aware ISO UTC timestamp")
    get.add_argument("--end", required=True, help="timezone-aware ISO UTC timestamp")
    get.add_argument("--exposure-id", type=int)
    get.add_argument("--limit", type=int)
    get.add_argument("--cache", type=Path, default=Path("cache"))
    get.add_argument("--output", type=Path, default=Path("out"))
    make = sub.add_parser("extract", help="measure calibrated, masked streak curves")
    make.add_argument("--manifest", type=Path, default=Path("out/manifest.json"))
    make.add_argument("--cache", type=Path, default=Path("cache"))
    make.add_argument("--output", type=Path, default=Path("out"))
    plot = sub.add_parser("plot", help="plot a tidy streak light curve")
    plot.add_argument("--curve", type=Path, required=True)
    plot.add_argument("--output", type=Path, required=True)
    plot.add_argument("--absolute-magnitude", action="store_true",
                      help="show header-calibrated magnitude instead of median-relative magnitude")
    plot.add_argument("--snr-cutoff", type=float, default=3.0,
                      help="minimum per-bin flux/uncertainty for filled markers and y limits (default: 3)")
    args = parser.parse_args(argv)
    try:
        if args.command == "fetch":
            result = fetch(args.norad_id, args.start, args.end, args.cache,
                           args.output, exposure_id=args.exposure_id, limit=args.limit)
        elif args.command == "extract":
            result = extract(args.manifest, args.cache, args.output)
        else:
            _plot(args.curve, args.output, absolute_magnitude=args.absolute_magnitude,
                  snr_cutoff=args.snr_cutoff)
            result = args.output
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"glint-geo: {exc}", file=sys.stderr)
        return 1
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
