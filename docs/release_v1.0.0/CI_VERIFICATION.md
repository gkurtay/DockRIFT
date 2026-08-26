# DockRIFT v1.0.0 CI verification

Release validation is performed by GitHub Actions on Python 3.10, 3.11, 3.12, and 3.13.

Each matrix job must complete:

- package installation from the public source tree;
- the 13-test regression and API suite; and
- `dockrift doctor` readiness checks.

After the matrix passes, a single build job creates the versioned wheel and source distribution directly from the tested public commit and records SHA256 digests for the release assets.

The tagged v1.0.0 release must be created only from a commit for which the complete matrix and build job are successful.
