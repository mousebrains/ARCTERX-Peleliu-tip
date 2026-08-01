#! /usr/bin/env python3

"""Loader and tidal analysis for the ARCTERX-2023 Peleliu Wake pressure array.

Twelve bottom-mounted wave gauges in 17-21 m of water spanning ~11 km: three at
the Peleliu tip (Pe1-3), five on Hydrographer Bank (HBN/S/E/W plus HBM on the
SeaSpider), four around Angaur (An1-4).  See PRESSURE_ANALYSIS.md for what the
array can and cannot measure; this module is only about reading it correctly.

Three things here exist because getting them wrong produced confident, wrong
answers -- and in two cases the wrong answer looked *better* by the obvious
metric:

**Declared time units are honoured, never assumed.**  Ten gauges say
"milliseconds since ..."; Pe2 (060465) correctly says "seconds since ...".
Hardcoding milliseconds turns Pe2's 31-day 1 Hz record into a phantom
45-minute record at 1000 Hz.  The physics catches it -- 2.7 m of tidal range
cannot happen in 45 minutes -- but only if you look.  ``xarray.open_dataset``
decodes CF time correctly and is immune; raw netCDF4 readers are not.

**Decimation block-averages, never strides.**  At 19 m an 8 s wave still
reaches the bottom at ~58 % amplitude, so 0.5 m of swell puts ~25 cm into
bottom pressure.  Striding to 1 min folds that into the tidal band: the
post-fit residual inflates from ~4 cm to 7.5-14 cm and the M2 phase spread
across the array becomes 149 deg, implying a 0.6 m/s tidal wave where shallow
water gives 13 m/s.

**Datasets are never closed.**  Pe2 segfaults inside HDF5 on
``Dataset.close()`` after a large read ("There are 1 HDF5 objects open!").
The reads themselves are correct, so every Dataset is kept referenced for the
life of the process and the OS reclaims the handles at exit.  Analysis scripts
are short-lived.  Do NOT reintroduce a context manager here.

The harmonic fit below is a thin local version.  Prefer
``wamos_tpw.tides.harmonic_fit``, which additionally refuses constituent sets
the record cannot resolve and reports the design-matrix condition number --
the only diagnostic that exposes a degenerate fit, because such a fit has a
*lower* residual than a correct one.

Usage::

    import pressure_array as pa
    t, dep, temp, lat, lon, rate = pa.load("Pe2", step_s=60.0)
    fit = pa.harmonic(t, dep)
    print(fit["amp"]["M2"], fit["var_explained"])
"""

from __future__ import annotations

import netCDF4, numpy as np, glob, os, re

from paths import ARCHIVE, DATA, have_archive

#: Eleven gauges live in "Pressure Sensors"; HBM sits on the SeaSpider mooring
#: in its own directory.  Search both so all twelve resolve by site name.
SEARCH_DIRS = (f"{ARCHIVE}/Pressure Sensors", f"{ARCHIVE}/Bank Seaspider")
BASE = SEARCH_DIRS[0]   # kept for backwards compatibility

#: The repository's vendored 1-minute block means.  Self-contained; used
#: whenever the requested step is >= 60 s, which covers every tidal, residual
#: and co-tidal analysis.  The full archive is needed only above the wave band.
VENDORED = f"{DATA}/pressure/1min"

#: Datasets are deliberately never closed -- see load().
_OPEN = []

# period in hours; the set resolvable in a ~26-31 day record
CONSTITUENTS = dict(
    M2=12.4206012, S2=12.0000000, N2=12.6583475,
    K1=23.9344697, O1=25.8193417,
    M4=6.2103006, MS4=6.1033393, M6=4.1402004,
    Mf=327.8599387, MSf=354.3670666,
)


def load(site, step_s=60.0, t0=None, t1=None, source="auto"):
    """One gauge block-averaged to ``step_s`` seconds.

    Args:
        site: Gauge name, e.g. ``"Pe1"``, ``"HBM"``.
        step_s: Output interval. Must be a multiple of 60 s when reading the
            vendored subset, which is already 1-minute means.
        t0, t1: Optional ``datetime64``-parseable bounds.
        source: ``"auto"`` prefers the vendored subset and falls back to the
            archive; ``"vendored"`` or ``"archive"`` force one.

    Returns:
        ``(times, dep, temp, lat, lon, rate_hz)``. ``rate_hz`` is the ORIGINAL
        instrument rate, not the decimated one.
    """
    use_vendored = (
        source == "vendored"
        or (source == "auto" and step_s >= 60.0
            and os.path.exists(f"{VENDORED}/{site}_1min.npz"))
    )
    if use_vendored:
        f = f"{VENDORED}/{site}_1min.npz"
        if not os.path.exists(f):
            raise FileNotFoundError(f"no vendored record for {site!r}: {f}")
        if step_s % 60 != 0:
            raise ValueError(
                f"vendored records are 1-minute means; step_s={step_s} is not a "
                "multiple of 60. Use source='archive' for finer output."
            )
        z = np.load(f, allow_pickle=False)
        t = z["time_ms"].astype("datetime64[ms]")
        dep = z["dep"].astype(float); temp = z["temp"].astype(float)
        k = int(step_s // 60)
        if k > 1:                       # further block-average, never stride
            n = (len(dep) // k) * k
            dep = dep[:n].reshape(-1, k).mean(axis=1)
            temp = temp[:n].reshape(-1, k).mean(axis=1)
            t = t[:n:k] + np.timedelta64(int(step_s * 500), "ms")
        sel = np.ones(len(t), bool)
        if t0 is not None: sel &= t >= np.datetime64(t0, "ms")
        if t1 is not None: sel &= t <= np.datetime64(t1, "ms")
        return (t[sel], dep[sel], temp[sel],
                float(z["lat"]), float(z["lon"]), float(z["source_rate_hz"]))

    if not have_archive():
        raise FileNotFoundError(
            f"{site!r} at step_s={step_s} needs the full archive, which is not "
            f"mounted at {ARCHIVE}.\n  Set ARCTERX_ARCHIVE, or use step_s>=60 "
            "to read the vendored 1-minute subset."
        )
    hits = [h for d in SEARCH_DIRS for h in glob.glob(f"{d}/*_{site}_clipped.nc")]
    if not hits:
        raise FileNotFoundError(
            f"no gauge named {site!r}; searched {list(SEARCH_DIRS)}"
        )
    f = hits[0]
    # Pe2 (060465) segfaults inside HDF5 on Dataset.close() after a large read
    # ("There are 1 HDF5 objects open!").  The reads themselves are correct, so
    # keep every Dataset referenced for the life of the process and never close
    # it.  Do NOT reintroduce a context manager here.
    ds = netCDF4.Dataset(f)
    _OPEN.append(ds)
    for v in ds.variables.values():
        v.set_auto_mask(False)
    tv = ds.variables["time"]
    # Honour the DECLARED unit.  Ten gauges say "milliseconds since ..."; Pe2
    # correctly says "seconds since ...".
    m = re.match(r"\s*(\w+)\s+since\s+(.+)", tv.units)
    if m is None:
        raise ValueError(f"cannot parse time units {tv.units!r}")
    unit, origin = m.group(1).lower(), m.group(2).strip()
    to_ms = {"seconds": 1000.0, "second": 1000.0, "milliseconds": 1.0,
             "microseconds": 1e-3, "minutes": 60000.0,
             "hours": 3600000.0, "days": 86400000.0}[unit]
    epoch = np.datetime64(origin.replace(" ", "T"), "ms")

    n_time = len(tv)
    head = np.asarray(tv[:1001], dtype=np.int64) * to_ms
    rate = 1000.0 / ((head[1000] - head[0]) / 1000.0)
    stride = max(int(round(step_s * rate)), 1)

    i0, i1 = 0, n_time
    if t0 is not None:
        i0 = max(int((np.datetime64(t0, "ms") - epoch)
                     / np.timedelta64(1, "ms") / 1000 * rate), 0)
    if t1 is not None:
        i1 = min(int((np.datetime64(t1, "ms") - epoch)
                     / np.timedelta64(1, "ms") / 1000 * rate), n_time)

    # BLOCK-AVERAGE, never stride -- see the module docstring.
    n = ((i1 - i0) // stride) * stride
    dep = np.asarray(ds.variables["dep"][i0:i0 + n], float
                     ).reshape(-1, stride).mean(axis=1)
    temp = np.asarray(ds.variables["temp"][i0:i0 + n], float
                      ).reshape(-1, stride).mean(axis=1)
    tsel = (np.asarray(tv[i0:i0 + n:stride], dtype=np.int64) * to_ms
            ).astype(np.int64) + int(stride / rate * 500)
    t = epoch + tsel.astype("timedelta64[ms]")
    lat = float(ds.variables["lat"][...])
    lon = float(ds.variables["lon"][...])
    return t, dep, temp, lat, lon, rate


PHASE_EPOCH = np.datetime64("2023-05-01T00:00:00", "ms")


def harmonic(t, y, constituents=CONSTITUENTS, trend=True, epoch=None):
    """Least-squares tidal fit. Returns dict of (amplitude m, phase deg) plus fit/residual."""
    ok = np.isfinite(y)
    # Phase referenced to a FIXED epoch shared by every sensor.  Referencing
    # to each record's own start makes the phases mutually meaningless.
    th = (t - (PHASE_EPOCH if epoch is None else epoch)) / np.timedelta64(1, "h")
    cols = [np.ones_like(th)]
    names = ["mean"]
    if trend:
        cols.append(th); names.append("trend")
    for nm, per in constituents.items():
        w = 2 * np.pi / per
        cols += [np.cos(w * th), np.sin(w * th)]
        names += [nm + "_c", nm + "_s"]
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A[ok], y[ok], rcond=None)
    fit = A @ coef
    res = y - fit
    out = {}
    for i, nm in enumerate(constituents):
        j = names.index(nm + "_c")
        c, s = coef[j], coef[j + 1]
        out[nm] = (float(np.hypot(c, s)), float(np.degrees(np.arctan2(-s, c)) % 360))
    return dict(amp=out, fit=fit, residual=res, coef=coef, names=names,
                var_explained=1 - np.nanvar(res[ok]) / np.nanvar(y[ok]))
