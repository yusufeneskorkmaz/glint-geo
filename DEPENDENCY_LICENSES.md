# Dependency license check for 0.1.2

The Python 3.12 environment resolved on 2026-09-27 was inspected with
`importlib.metadata` before the first push. The direct runtime dependencies
reported: Astropy BSD-3-Clause, HTTPX BSD, Matplotlib PSF, NumPy BSD-3-Clause
with 0BSD/MIT/Zlib/CC0 components, SciPy BSD, and Skyfield MIT. Skyfield's
SGP4 dependency reported MIT. The test-only dependency pytest reported MIT.

The installed transitive dependency metadata showed MIT, BSD, Apache, PSF,
MPL-2.0, or similarly permissive licenses; no GPL license was reported. The
two packages whose installed metadata did not state an identifier were checked
against their publishers' sources: [fontTools is MIT](https://github.com/fonttools/fonttools)
and [jplephem is MIT](https://github.com/brandon-rhodes/python-jplephem/blob/master/setup.py).
The repository vendors no dependency source or binaries.
