# Building DockRIFT release artifacts

The reference v1.0.0 wheel was built with Python 3.13 and `SOURCE_DATE_EPOCH=1787702400` (2026-08-26 00:00:00 UTC).

```bash
python -m pip install --upgrade pip setuptools
python -m pytest -q
rm -rf build dist src/dockrift.egg-info
mkdir -p dist
SOURCE_DATE_EPOCH=1787702400 python -m pip wheel . --no-deps --no-build-isolation -w dist
SOURCE_DATE_EPOCH=1787702400 python - <<'PY'
from setuptools.build_meta import build_sdist
build_sdist('dist')
PY
sha256sum dist/* | sed 's#dist/##' > RELEASE_SHA256SUMS_v1.0.0.txt
```

The release wheel should then be installed into an isolated target or environment and the supplied tests, `dockrift doctor`, and VDR API smoke workflow rerun before publication.
