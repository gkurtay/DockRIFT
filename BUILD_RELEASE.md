# Building DockRIFT release artifacts

The canonical DockRIFT v1.0.0 distributions are built by GitHub Actions on Python 3.13 with `SOURCE_DATE_EPOCH=1787702400` (2026-08-26 00:00:00 UTC), `pip==26.2.1`, and `setuptools==84.0.0`.

The wheel is built directly from the tested public source tree. The setuptools source distribution is then normalized to a deterministic tar/gzip representation by fixing file order, mtimes, ownership metadata, tar format, and gzip header timestamps.

```bash
python -m pip install "pip==26.2.1" "setuptools==84.0.0"
rm -rf build dist src/dockrift.egg-info .sdist-normalize
mkdir -p dist

SOURCE_DATE_EPOCH=1787702400 \
  python -m pip wheel . --no-deps --no-build-isolation -w dist

SOURCE_DATE_EPOCH=1787702400 python - <<'PY'
from setuptools.build_meta import build_sdist
build_sdist('dist')
PY

mkdir -p .sdist-normalize
tar -xzf dist/dockrift-1.0.0.tar.gz -C .sdist-normalize
rm dist/dockrift-1.0.0.tar.gz
(
  cd .sdist-normalize
  LC_ALL=C tar \
    --sort=name \
    --mtime='@1787702400' \
    --owner=0 \
    --group=0 \
    --numeric-owner \
    --format=gnu \
    -cf - dockrift-1.0.0
) | gzip -n > dist/dockrift-1.0.0.tar.gz
rm -rf .sdist-normalize

sha256sum dist/* | sed 's#dist/##'
```

The expected v1.0.0 hashes are recorded in `RELEASE_SHA256SUMS_v1.0.0.txt`. Before publication, the complete Python 3.10–3.13 test matrix, `dockrift doctor`, and the release build must pass from the exact tagged commit. The distributions attached to the GitHub release must match the recorded SHA256 values.
