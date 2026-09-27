"""Command-line interface for acquisition, extraction, and plotting."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys

from .pipeline import extract, fetch


def _plot(curve: Path, output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = list(csv.DictReader(curve.open()))
    good = [row for row in rows if row["quality"] == "OK" and row["ZTF_mag_no_color"]]
    if not good:
        raise ValueError("curve has no calibrated, unmasked samples")
    t = [float(row["seconds_from_shutter_open"]) for row in good]
    mag = [float(row["ZTF_mag_no_color"]) for row in good]
    error = [float(row["mag_uncertainty"]) for row in good]
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.errorbar(t, mag, yerr=error, fmt=".", color="navy", ecolor="steelblue", capsize=0)
    ax.invert_yaxis()
    ax.set(xlabel="Seconds from shutter open", ylabel="ZTF magnitude (no color correction)")
    ax.grid(alpha=0.2)
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
    args = parser.parse_args(argv)
    try:
        if args.command == "fetch":
            result = fetch(args.norad_id, args.start, args.end, args.cache,
                           args.output, exposure_id=args.exposure_id, limit=args.limit)
        elif args.command == "extract":
            result = extract(args.manifest, args.cache, args.output)
        else:
            _plot(args.curve, args.output)
            result = args.output
    except (ValueError, RuntimeError, OSError) as exc:
        print(f"glint-geo: {exc}", file=sys.stderr)
        return 1
    print(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
