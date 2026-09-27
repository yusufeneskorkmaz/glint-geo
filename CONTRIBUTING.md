# Contributing

Open an issue before a substantial change. Use Python 3.11 or newer and run
`python -m pytest` locally. Add synthetic offline tests for changes to geometry,
masking, or calibration. Mark any network test `@pytest.mark.live`.

Never commit credentials, orbital elements, TLEs, downloaded images, field-star
catalogs, or generated light curves. Keep acknowledgments in `DATA_POLICY.md`
when publishing results. Contributions use the BSD-3-Clause license.
