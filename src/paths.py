#! /usr/bin/env python3

"""Where the data lives.

Two roots, and the distinction matters:

``DATA``
    The repository's vendored subset (``data/``). Self-contained: every
    analysis in PRESSURE_ANALYSIS.md and DRIFTER_ANALYSIS.md runs from it with
    no external volume mounted. Pressure records here are 1-minute block means
    (see ``tools/make_data_subset.py``), which is everything the tidal,
    residual and co-tidal work needs.

``ARCHIVE``
    The full SeaChest archive, ~2.6 GB, needed only for work above the wave
    band at gauges other than the two vendored full-rate segments. Absent on
    most machines; code must degrade gracefully rather than crash.

Override either with an environment variable::

    export ARCTERX_DATA=/somewhere/else/data
    export ARCTERX_ARCHIVE=/Volumes/SeaChest/ARCTERX/2023/Wake
"""

from __future__ import annotations

import os

__all__ = ["REPO", "DATA", "ARCHIVE", "have_archive", "drifter_nc", "require"]

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

DATA = os.path.abspath(os.environ.get("ARCTERX_DATA", os.path.join(REPO, "data")))

ARCHIVE = os.environ.get("ARCTERX_ARCHIVE",
                         "/Volumes/SeaChest/ARCTERX/2023/Wake")

#: The four drifters that followed the 2023-05-22 eddy.
DRIFTERS = ("mwb458d02", "mwb788d01", "mwb790d01", "mwb793d02")


def have_archive() -> bool:
    """True when the full SeaChest archive is reachable."""
    return os.path.isdir(ARCHIVE)


def drifter_nc(name: str) -> str:
    """Path to one regenerated drifter timeseries."""
    return os.path.join(DATA, "drifters", f"{name}_gps_timeseries.nc")


def require(path: str, what: str = "") -> str:
    """Return ``path``, or raise with a message that says how to obtain it."""
    if os.path.exists(path):
        return path
    hint = (f"\n  It lives in the full archive ({ARCHIVE}), which is not mounted."
            "\n  Either mount it, set ARCTERX_ARCHIVE, or rebuild the vendored"
            "\n  subset with tools/make_data_subset.py."
            if not have_archive() else "")
    raise FileNotFoundError(f"missing {what or 'data'}: {path}{hint}")
