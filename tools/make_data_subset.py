#! /usr/bin/env python3

"""Build the repository's vendored ``data/`` subset from the SeaChest archive.

The full archive is ~2.6 GB, almost all of it bottom-pressure records sampled
at 2-16 Hz for a month. Every tidal, residual, co-tidal and gradient result in
PRESSURE_ANALYSIS.md comes from **1-minute block means**, so the vendored copy
is decimated by a factor of ~180 with nothing lost above the wave band. Two
short high-rate segments are kept so the aliasing and cosh(kh) lessons stay
runnable without the archive.

This script exists so the trimming is auditable and repeatable rather than a
one-off nobody can reproduce. Run it only when regenerating ``data/``; the
repository ships the output.

    python3 tools/make_data_subset.py --archive /Volumes/SeaChest/ARCTERX/2023/Wake

Decimation is by BLOCK AVERAGE, never striding -- see PRESSURE_ANALYSIS.md
section 3.2 for what striding costs.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: All twelve gauges: eleven in "Pressure Sensors", HBM on the SeaSpider.
GAUGES = ("Pe1", "Pe2", "Pe3", "HBN", "HBS", "HBE", "HBW", "HBM",
          "An1", "An2", "An3", "An4")

#: Short full-rate windows kept for the teaching demos.
HIGHRATE = (
    ("Pe3", "2023-05-22T00:00", "2023-05-22T06:00", "16 Hz, cosh(kh) wave demo"),
    ("HBN", "2023-05-22T00:00", "2023-05-23T00:00", "2 Hz, stride-vs-average demo"),
)

#: Bounding box covering the drifter tracks and every sensor, with margin.
BBOX = dict(lat0=6.88, lat1=7.05, lon0=134.05, lon1=134.28)


def _open(path):
    import netCDF4
    ds = netCDF4.Dataset(path)
    for v in ds.variables.values():
        v.set_auto_mask(False)
    return ds


def _gauge_path(archive: str, site: str) -> str:
    import glob
    for d in (f"{archive}/Pressure Sensors", f"{archive}/Bank Seaspider"):
        hits = glob.glob(f"{d}/*_{site}_clipped.nc")
        if hits:
            return hits[0]
    raise FileNotFoundError(f"no gauge {site!r} under {archive}")


def _time_scale(units: str) -> tuple[float, np.datetime64]:
    """Honour the DECLARED unit -- Pe2 says 'seconds', the rest 'milliseconds'."""
    m = re.match(r"\s*(\w+)\s+since\s+(.+)", units)
    if m is None:
        raise ValueError(f"cannot parse time units {units!r}")
    to_ms = {"seconds": 1000.0, "milliseconds": 1.0, "microseconds": 1e-3,
             "minutes": 60000.0, "hours": 3600000.0, "days": 86400000.0}
    return to_ms[m.group(1).lower()], np.datetime64(m.group(2).strip().replace(" ", "T"), "ms")


def pressure_1min(archive: str, out_dir: str) -> None:
    """All twelve gauges, block-averaged to 1 minute."""
    os.makedirs(out_dir, exist_ok=True)
    for site in GAUGES:
        path = _gauge_path(archive, site)
        ds = _open(path)
        tv = ds.variables["time"]
        to_ms, epoch = _time_scale(tv.units)
        head = np.asarray(tv[:1001], dtype=np.int64) * to_ms
        rate = 1000.0 / ((head[1000] - head[0]) / 1000.0)
        stride = max(int(round(60.0 * rate)), 1)
        n = (len(tv) // stride) * stride
        dep = np.asarray(ds.variables["dep"][:n], float).reshape(-1, stride).mean(axis=1)
        tmp = np.asarray(ds.variables["temp"][:n], float).reshape(-1, stride).mean(axis=1)
        tms = ((np.asarray(tv[0:n:stride], dtype=np.int64) * to_ms).astype(np.int64)
               + int(stride / rate * 500))
        np.savez_compressed(
            os.path.join(out_dir, f"{site}_1min.npz"),
            time_ms=(epoch.astype(np.int64) + tms).astype(np.int64),
            dep=dep.astype(np.float32), temp=tmp.astype(np.float32),
            lat=float(ds.variables["lat"][...]), lon=float(ds.variables["lon"][...]),
            source_rate_hz=rate, source_file=os.path.basename(path),
        )
        print(f"  {site:4s} {rate:5.1f} Hz -> {len(dep):6d} one-minute means")
        # NB: never ds.close() -- Pe2 segfaults inside HDF5 on close after a
        # large read.  See PRESSURE_ANALYSIS.md section 3.6.


def pressure_highrate(archive: str, out_dir: str) -> None:
    """Short full-rate windows for the aliasing and wave-attenuation demos."""
    os.makedirs(out_dir, exist_ok=True)
    for site, t0, t1, why in HIGHRATE:
        ds = _open(_gauge_path(archive, site))
        tv = ds.variables["time"]
        to_ms, epoch = _time_scale(tv.units)
        head = np.asarray(tv[:1001], dtype=np.int64) * to_ms
        rate = 1000.0 / ((head[1000] - head[0]) / 1000.0)
        i0 = int((np.datetime64(t0, "ms") - epoch) / np.timedelta64(1, "ms") / 1000 * rate)
        i1 = int((np.datetime64(t1, "ms") - epoch) / np.timedelta64(1, "ms") / 1000 * rate)
        i0, i1 = max(i0, 0), min(i1, len(tv))
        np.savez_compressed(
            os.path.join(out_dir, f"{site}_{t0[:10]}_fullrate.npz"),
            time_ms=(epoch.astype(np.int64)
                     + (np.asarray(tv[i0:i1], dtype=np.int64) * to_ms).astype(np.int64)),
            dep=np.asarray(ds.variables["dep"][i0:i1], np.float32),
            lat=float(ds.variables["lat"][...]), lon=float(ds.variables["lon"][...]),
            rate_hz=rate, purpose=why,
        )
        print(f"  {site:4s} {rate:5.1f} Hz x {(i1-i0)/rate/3600:.1f} h -> {i1-i0} samples ({why})")


def adcp_c05(archive: str, out_path: str) -> None:
    """C05 depth-averaged current and depth, WITHOUT the rotation applied.

    The rotation is left to the analysis so the raw instrument frame stays
    visible; see PRESSURE_ANALYSIS.md section 4 for the correction.
    """
    import scipy.io as sio
    a = sio.loadmat(f"{archive}/Bank ADCP/C05_2023_proc.mat",
                    squeeze_me=True, struct_as_record=False)["adcpr"]
    t = ((a.mtime - 719529.0) * 86400 * 1000).astype(np.int64)
    E = np.asarray(a.east_vel, float)
    N = np.asarray(a.north_vel, float)
    with np.errstate(invalid="ignore"):
        Em, Nm = np.nanmean(E, axis=0), np.nanmean(N, axis=0)
    dep = np.asarray(a.depth, float)
    ok = np.isfinite(Em) & np.isfinite(Nm) & (dep > 5)   # 0 m = out of water
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    np.savez_compressed(
        out_path, time_ms=t[ok],
        east=Em[ok].astype(np.float32), north=Nm[ok].astype(np.float32),
        depth=dep[ok].astype(np.float32),
        temperature=np.asarray(a.temperature, np.float32)[ok],
        ranges=np.asarray(a.ranges, np.float32),
        lat=6.9302, lon=134.1994,
        note="depth-averaged; rotation NOT applied; 12-min ensembles",
    )
    print(f"  C05 {ok.sum()} of {len(t)} ensembles kept (depth>5 m mask)")


def bathymetry(archive: str, out_path: str) -> None:
    """Crop the 25 m DEM to the working box, keep Z only, store float32."""
    import netCDF4
    ds = _open(f"{archive}/bathy/Angaur_Peleliu_25m.nc")
    lat = np.asarray(ds.variables["lat"][:], float)
    lon = np.asarray(ds.variables["lon"][:], float)
    Z = np.asarray(ds.variables["Z"][:], float)
    inbox = ((lat >= BBOX["lat0"]) & (lat <= BBOX["lat1"])
             & (lon >= BBOX["lon0"]) & (lon <= BBOX["lon1"]))
    rows = np.where(inbox.any(axis=1))[0]
    cols = np.where(inbox.any(axis=0))[0]
    r0, r1, c0, c1 = rows[0], rows[-1] + 1, cols[0], cols[-1] + 1
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with netCDF4.Dataset(out_path, "w", format="NETCDF4") as out:
        out.createDimension("y", r1 - r0)
        out.createDimension("x", c1 - c0)
        for name, arr, units, long_name in (
                ("lat", lat[r0:r1, c0:c1], "degrees_north", "latitude"),
                ("lon", lon[r0:r1, c0:c1], "degrees_east", "longitude"),
                ("Z", Z[r0:r1, c0:c1], "m", "elevation, positive up")):
            v = out.createVariable(name, "f4", ("y", "x"), zlib=True, complevel=5)
            v[:] = arr.astype(np.float32)
            v.units = units
            v.long_name = long_name
        out.title = "Angaur/Peleliu 25 m DEM, cropped"
        out.source = "Angaur_Peleliu_25m.nc, cropped and cast to float32"
        out.comment = ("DEM is unreliable over steep reef: within 150 m of the "
                       "Peleliu and Angaur gauges the relief is 18-59 m. See "
                       "PRESSURE_ANALYSIS.md section 6.")
    print(f"  bathymetry cropped to {r1-r0} x {c1-c0} cells")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--archive", default="/Volumes/SeaChest/ARCTERX/2023/Wake")
    p.add_argument("--out", default=os.path.join(REPO, "data"))
    p.add_argument("--skip", nargs="*", default=[],
                   choices=["pressure", "highrate", "adcp", "bathy"])
    a = p.parse_args(argv)
    if not os.path.isdir(a.archive):
        print(f"archive not found: {a.archive}", file=sys.stderr)
        return 1
    if "pressure" not in a.skip:
        print("pressure gauges -> 1-minute block means")
        pressure_1min(a.archive, os.path.join(a.out, "pressure", "1min"))
    if "highrate" not in a.skip:
        print("full-rate demo segments")
        pressure_highrate(a.archive, os.path.join(a.out, "pressure", "highrate"))
    if "adcp" not in a.skip:
        print("C05 ADCP")
        adcp_c05(a.archive, os.path.join(a.out, "adcp", "c05_depth_avg.npz"))
    if "bathy" not in a.skip:
        print("bathymetry")
        bathymetry(a.archive, os.path.join(a.out, "bathy", "angaur_peleliu_crop.nc"))
    print("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
