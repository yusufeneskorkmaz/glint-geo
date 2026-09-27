# glint-geo

`glint-geo` extracts calibrated, sub-second light curves from public Zwicky
Transient Facility science-image streaks at predicted locations of catalogued
geosynchronous objects. It uses SGP4 predictions and the public IRSA metadata,
science-cutout, mask-cutout, and Gaia-catalog services. It is a photometry tool,
not an attitude-state estimator.

## Installation

Python 3.11 or newer is required.

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install -e .
```

## Five-minute quickstart

The GOES 12 exposure below is a public exposure ID identified by the GEO
survey. Historical dates require your own Space-Track account; set the two
environment variables in your shell or credential manager. Do not put them in
a notebook or `.env` file.

```bash
glint-geo fetch --norad-id 26871 \
  --start 2019-06-15T04:58:00Z --end 2019-06-15T05:00:00Z \
  --exposure-id 89520732 --limit 1
glint-geo extract
glint-geo plot --curve out/26871_89520732_c05q2.csv \
  --output out/goes12.png
```

`fetch` writes a provenance manifest and downloads public cutouts to `cache/`.
`extract` writes one tidy CSV per streak plus calibration, URLs, hashes, and
quality counts in `out/manifest.json`. Each CSV contains flux, uncertainty,
color-uncorrected ZTF magnitude, time from shutter open, and a per-bin quality
flag. Cache and output folders are excluded from Git. Use `--cache` and
`--output` to choose other locations. Run `glint-geo --help` for the CLI.

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

## Tests and citation

`python -m pytest` runs synthetic offline tests. Live tests are marked `live`
and skipped by default. See `examples/goes12_quickstart.ipynb` for the same
public-data workflow in a notebook.

Use the citation in `CITATION.cff`. For ZTF and IRSA acknowledgment text and
orbital-data handling, see `DATA_POLICY.md`. The public ZTF image service is
identified by [DOI 10.26131/IRSA539](https://doi.org/10.26131/IRSA539).
