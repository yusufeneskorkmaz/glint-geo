# glint-geo

`glint-geo` extracts calibrated, sub-second light curves from public Zwicky
Transient Facility science-image streaks at predicted locations of catalogued
geosynchronous objects. It uses SGP4 predictions and the public IRSA metadata,
science-cutout, mask-cutout, and Gaia-catalog services. It is a photometry tool,
not an attitude-state estimator.

## Installation

Python 3.11 or newer is required. Install the test extra before running tests:

```bash
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -e ".[test]"
python3 -m pytest
```

## Five-minute quickstart

The Figure 2 commands below fetch one public GOES 12 exposure, extract its
streak light curve, and plot it. Historical dates require your own Space-Track
account. Set `SPACETRACK_USER` and `SPACETRACK_PASSWORD` in the shell
environment before running the commands. Never write them to files, including
notebooks or `.env` files.

## Reproduce Figure 2 of the technical note

Figure 2 and its interpretation are in the [included technical note](docs/technical_note.pdf).
Exposure 115845983 was acquired on the 2020-03-03 Palomar night, at
2020-03-04 11:02 UTC. Run these commands from the repository root:

```bash
python3 -m glint_geo.cli fetch --norad-id 26871 \
  --start 2020-03-04T11:02:00Z --end 2020-03-04T11:03:00Z \
  --exposure-id 115845983 --limit 1
python3 -m glint_geo.cli extract
python3 -m glint_geo.cli plot --curve out/26871_115845983_c07q2.csv \
  --output out/goes12_figure2.png
```

`fetch` writes a provenance manifest and downloads public cutouts to `cache/`.
`extract` writes one tidy CSV per streak plus calibration, URLs, hashes, and
quality counts in `out/manifest.json`. Each CSV contains flux, uncertainty,
color-uncorrected ZTF magnitude, time from shutter open, and a per-bin quality
flag. Cache and output folders are excluded from Git. Use `--cache` and
`--output` to choose other locations. Run `python3 -m glint_geo.cli --help`
for the CLI.

The fallback CelesTrak query supplies *current* elements only. It is available
for recent observation epochs when Space-Track variables are absent. Historical
queries fail clearly without Space-Track access. An element more than three
days from an observation is never used. Cutout matching requires both predicted
trail endpoints to lie within one science quadrant.

## Interpretation

Resolved streaks can measure repeated brightness structure from sufficiently
fast rotators within a single exposure. Sparse visits across nights cannot
reliably link the phase of slow tumblers or changing spin states. Header
zero points are field-star calibrations; reported magnitudes have no satellite
color correction, and the header RMS is not the full photometric uncertainty.
Mask and bright-star crossings are excluded by segment, with quality flags for
each bin. Inspect the cutout and mask before treating a curve as physical.

`plot` defaults to relative magnitude: it subtracts the median magnitude of
unmasked bins with per-bin flux/uncertainty **SNR ≥ 3**, matching the relative
scale of the technical note figures. Filled points meet that cutoff; lower-SNR
bins use hollow markers and do not set the y-axis range. Hollow markers at an
axis edge represent off-scale or non-positive-flux bins. Shaded intervals
identify bad-pixel and bright-star masks in the legend. Use
`--absolute-magnitude` for the header-calibrated, color-uncorrected scale or
`--snr-cutoff` to change the display cutoff. These options affect plots only;
they do not change the extracted CSV or provenance manifest.

## Tests and citation

`python3 -m pytest` runs synthetic offline tests. Live tests are marked `live`
and skipped by default. See `examples/goes12_quickstart.ipynb` for the same
public-data workflow in a notebook.

Technical note citation: Korkmaz, Y. E. (2026). *Archival ZTF streak photometry
of retired GEO satellites.* Technical note, [docs/technical_note.pdf](docs/technical_note.pdf).
The PDF states a CC-BY-4.0 license; the package code is BSD-3-Clause.
For the software, use `CITATION.cff`. For ZTF and IRSA acknowledgment text and
orbital-data handling, see `DATA_POLICY.md`. The public ZTF image service is
identified by [DOI 10.26131/IRSA539](https://doi.org/10.26131/IRSA539).
