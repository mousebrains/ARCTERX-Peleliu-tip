#! /usr/bin/env python3

"""Reader for CORDC Miniature Wave Buoy (mwb) raw GPS ``.dat`` files.

FILE FORMAT
-----------
Each ``MWB_<yyyymmddHHMMSS>.dat`` is one sampling burst ("dwell"): exactly
2048 fixed-size records of 85 bytes, little-endian, no header and no footer
(174080 bytes total).  The filename stamp is the burst *command* time, which
runs 5.0-6.0 s ahead of the first record; use the record fields for time.

The 85-byte record is a truncated UBX-NAV-PVT: field order and byte offsets
follow u-blox NAV-PVT exactly, but the fields that u-blox stores as scaled
integers (lon, lat, height, hMSL, velN, velE, velD, gSpeed, headMot) are
stored here as IEEE-754 float32 already in physical units.  Offsets 78-84 are
zero padding.

  off  size  type   name      units as stored
    0     4  u4     iTOW      GPS time of week, ms
    4     2  u2     year      UTC
    6     1  u1     month
    7     1  u1     day
    8     1  u1     hour
    9     1  u1     min
   10     1  u1     sec
   11     1  u1     valid     bitfield
   12     4  u4     tAcc      time accuracy estimate, ns
   16     4  i4     nano      sub-second residual of UTC time, ns (may be < 0)
   20     1  u1     fixType   3 = 3D fix
   21     1  u1     flags     bit0 gnssFixOK, bit1 diffSoln
   22     1  u1     flags2
   23     1  u1     numSV     satellites used in solution
   24     4  f4     lon       degrees east
   28     4  f4     lat       degrees north
   32     4  f4     height    height above WGS-84 ellipsoid, m
   36     4  f4     hMSL      height above mean sea level, m
   40     4  u4     hAcc      horizontal accuracy estimate, mm
   44     4  u4     vAcc      vertical accuracy estimate, mm
   48     4  f4     velN      north velocity, m/s
   52     4  f4     velE      east velocity, m/s
   56     4  f4     velD      DOWN velocity, m/s
   60     4  f4     gSpeed    2-D ground speed, m/s
   64     4  f4     headMot   heading of motion, degrees
   68     4  u4     sAcc      speed accuracy estimate, mm/s
   72     4  u4     headAcc   heading accuracy estimate, 1e-5 degrees
   76     2  u2     pDOP      position dilution of precision, 0.01
   78     7  --     (zero padding)

PROVENANCE OF THE FIELD IDENTIFICATIONS
---------------------------------------
Confirmed bit-for-bit against ``mwb458d02_gps_timeseries.nc`` over all 84
delivered dwells (172032 samples), max |difference| exactly 0:
    time (calendar fields + nano), lat, lon, u == velE, v == velN,
    w == -velD, numsats == numSV.
Confirmed to the netCDF's own quantization:
    sog == round(gSpeed*10)/10  (100.0000% of samples)
    cog == round(headMot*10)/10 ( 99.9948%; 9 samples differ by one
                                  0.1 deg count, float32 tie-breaking)
    pdop: the netCDF stores round(pDOP/100) as a byte (100% match), i.e.
          it keeps only integer PDOP.  The raw .dat keeps 0.01 resolution.

Inferred from the UBX-NAV-PVT convention rather than from ground truth --
the netCDF does not carry these -- but supported by independent checks:
    iTOW, nano, tAcc, valid, fixType, flags, flags2, height, hMSL,
    hAcc, vAcc, sAcc, headAcc.
  * iTOW - (UTC + 18 s leap seconds), reduced modulo the GPS week, is
    identically zero across all 116 non-empty files.
  * height - hMSL = 64.35 m, matching the EGM96 geoid undulation at Palau
    (~64 m), so the two height fields and their order are right.
  * nano makes the sample cadence exactly 500 ms.
  * bytes 78-84 are zero in all 237568 records, which also pins the record
    length: any misalignment would scramble that column.
The mm / mm/s / 1e-5 deg scalings for hAcc, vAcc, sAcc and headAcc are the
u-blox convention and give physically sensible magnitudes (1.0-1.3 m,
1.4-1.7 m, 0.11-0.19 m/s, 4.8-9.8 deg), but are NOT independently verified.

KNOWN LIMITATIONS OF THE SOURCE DATA
------------------------------------
* lon/lat are float32, so their resolution depends on magnitude.  At this
  deployment (7 N, 134 E) one float32 ulp is 1.53e-5 deg in longitude
  (1.69 m) but only 4.8e-7 deg in latitude (0.05 m) -- longitude is 32x
  coarser, purely because 134 and 7 sit in different binades.  Measured
  spacing between realized values in the ARCTERX-2023 data is 1.53e-5 deg
  (lon) and 9.5e-7 deg (lat).  This is a property of the file format, not
  of this reader; the delivered netCDF inherits it too.
* headMot is the receiver's filtered course.  It tracks
  atan2(velE, velN) to a median 0.76 deg / p90 1.9 deg, but diverges
  arbitrarily when gSpeed is near zero, where course is ill-defined.
* Receiver clock jitter of +/-1 ms appears occasionally, so consecutive
  sample spacing is 499-501 ms rather than exactly 500 ms.

USAGE
-----
    import mwb_dat

    rec = mwb_dat.read_file("raw_data/458/20230522_00/MWB_20230522000422.dat")
    ts  = mwb_dat.read_tree("raw_data/458")     # flat, sorted, de-duplicated
    print(ts["time"][0], ts["lat"][0], ts["u"][0])

Command line:
    python3 mwb_dat.py raw_data/458                    # summary
    python3 mwb_dat.py raw_data/458 -o out.npz         # numpy archive
    python3 mwb_dat.py raw_data/458 -o out.csv         # csv
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
import warnings

import numpy as np

__all__ = [
    "RECORD_DTYPE", "RECORD_SIZE", "SAMPLES_PER_DWELL",
    "read_file", "read_raw", "read_tree", "decode", "utc_ms", "save",
]

RECORD_SIZE = 85
SAMPLES_PER_DWELL = 2048

#: Byte-for-byte layout of one record, as stored on disk.
RECORD_DTYPE = np.dtype([
    ("iTOW",    "<u4"),
    ("year",    "<u2"),
    ("month",   "u1"),
    ("day",     "u1"),
    ("hour",    "u1"),
    ("min",     "u1"),
    ("sec",     "u1"),
    ("valid",   "u1"),
    ("tAcc",    "<u4"),
    ("nano",    "<i4"),
    ("fixType", "u1"),
    ("flags",   "u1"),
    ("flags2",  "u1"),
    ("numSV",   "u1"),
    ("lon",     "<f4"),
    ("lat",     "<f4"),
    ("height",  "<f4"),
    ("hMSL",    "<f4"),
    ("hAcc",    "<u4"),
    ("vAcc",    "<u4"),
    ("velN",    "<f4"),
    ("velE",    "<f4"),
    ("velD",    "<f4"),
    ("gSpeed",  "<f4"),
    ("headMot", "<f4"),
    ("sAcc",    "<u4"),
    ("headAcc", "<u4"),
    ("pDOP",    "<u2"),
    ("_pad",    "V7"),
])
assert RECORD_DTYPE.itemsize == RECORD_SIZE

GPS_UTC_LEAP_SECONDS = 18      # valid from 2017-01-01 until the next leap second
_MS_PER_WEEK = 7 * 86_400_000


# ----------------------------------------------------------------- time ----

def _days_from_civil(y, m, d):
    """Days since 1970-01-01 from proleptic-Gregorian y/m/d (Hinnant's algorithm)."""
    y = np.asarray(y, dtype=np.int64) - (np.asarray(m, dtype=np.int64) <= 2)
    m = np.asarray(m, dtype=np.int64)
    d = np.asarray(d, dtype=np.int64)
    era = np.where(y >= 0, y, y - 399) // 400
    yoe = y - era * 400
    doy = (153 * (m + np.where(m > 2, -3, 9)) + 2) // 5 + d - 1
    doe = yoe * 365 + yoe // 4 - yoe // 100 + doy
    return era * 146_097 + doe - 719_468


def utc_ms(raw):
    """UTC milliseconds since 1970-01-01 for each record.

    This is the receiver's UTC calendar time plus the ``nano`` residual,
    rounded to the millisecond -- the exact quantity the CORDC netCDF
    ``time`` variable carries.
    """
    days = _days_from_civil(raw["year"], raw["month"], raw["day"])
    ms = (days * 86_400_000
          + raw["hour"].astype(np.int64) * 3_600_000
          + raw["min"].astype(np.int64) * 60_000
          + raw["sec"].astype(np.int64) * 1_000)
    return ms + np.round(raw["nano"] / 1e6).astype(np.int64)


# ------------------------------------------------------------- decoding ----

#: Layout of the decoded output.  Fields that exist in the delivered netCDF
#: keep their netCDF names (time, lat, lon, u, v, w, sog, cog, numsats, pdop).
OUTPUT_DTYPE = np.dtype([
    ("time",      "<i8"),   # ms since 1970-01-01 UTC
    ("datetime",  "<M8[ms]"),
    ("lat",       "<f4"),   # degrees north  (float32 as stored)
    ("lon",       "<f4"),   # degrees east   (float32 as stored)
    ("u",         "<f4"),   # eastward velocity, m/s      (= velE)
    ("v",         "<f4"),   # northward velocity, m/s     (= velN)
    ("w",         "<f4"),   # upward velocity, m/s        (= -velD)
    ("sog",       "<f4"),   # speed over ground, m/s      (= gSpeed)
    ("cog",       "<f4"),   # course over ground, degrees (= headMot)
    ("numsats",   "u1"),
    ("pdop",      "<f4"),   # dimensionless, 0.01 resolution
    ("height",    "<f4"),   # above WGS-84 ellipsoid, m
    ("hmsl",      "<f4"),   # above mean sea level, m
    ("h_acc",     "<f4"),   # horizontal accuracy estimate, m
    ("v_acc",     "<f4"),   # vertical accuracy estimate, m
    ("s_acc",     "<f4"),   # speed accuracy estimate, m/s
    ("head_acc",  "<f4"),   # heading accuracy estimate, degrees
    ("t_acc",     "<u4"),   # time accuracy estimate, ns
    ("nano",      "<i4"),   # sub-second residual of UTC time, ns
    ("itow",      "<u4"),   # GPS time of week, ms
    ("fix_type",  "u1"),    # 3 = 3D fix
    ("flags",     "u1"),
    ("flags2",    "u1"),
    ("valid",     "u1"),
    ("dwell",     "<i4"),   # index of the source burst, in file order
])


def decode(raw, dwell=-1):
    """Convert raw records to physical units.

    Parameters
    ----------
    raw : structured ndarray of ``RECORD_DTYPE``
    dwell : int
        Value to store in the ``dwell`` column.

    Returns
    -------
    structured ndarray of ``OUTPUT_DTYPE``
    """
    out = np.empty(len(raw), dtype=OUTPUT_DTYPE)
    t = utc_ms(raw)
    out["time"] = t
    out["datetime"] = t.astype("<M8[ms]")
    out["lat"] = raw["lat"]
    out["lon"] = raw["lon"]
    out["u"] = raw["velE"]
    out["v"] = raw["velN"]
    out["w"] = -raw["velD"]          # velD is positive DOWN; w is positive up
    out["sog"] = raw["gSpeed"]
    out["cog"] = raw["headMot"]
    out["numsats"] = raw["numSV"]
    out["pdop"] = raw["pDOP"] * 0.01
    out["height"] = raw["height"]
    out["hmsl"] = raw["hMSL"]
    out["h_acc"] = raw["hAcc"] * 1e-3        # mm -> m
    out["v_acc"] = raw["vAcc"] * 1e-3        # mm -> m
    out["s_acc"] = raw["sAcc"] * 1e-3        # mm/s -> m/s
    out["head_acc"] = raw["headAcc"] * 1e-5  # 1e-5 deg -> deg
    out["t_acc"] = raw["tAcc"]
    out["nano"] = raw["nano"]
    out["itow"] = raw["iTOW"]
    out["fix_type"] = raw["fixType"]
    out["flags"] = raw["flags"]
    out["flags2"] = raw["flags2"]
    out["valid"] = raw["valid"]
    out["dwell"] = dwell
    return out


# --------------------------------------------------------------- readers ---

def read_raw(path, strict=False):
    """Read one ``.dat`` and return the undecoded records (``RECORD_DTYPE``).

    Returns an empty array for a zero-length file.  A file whose length is
    not a whole number of records is truncated to the last complete record
    (with a warning) unless ``strict``.
    """
    size = os.path.getsize(path)
    if size == 0:
        if strict:
            raise ValueError(f"{path}: empty file")
        warnings.warn(f"{path}: empty file, skipped", stacklevel=2)
        return np.empty(0, dtype=RECORD_DTYPE)
    extra = size % RECORD_SIZE
    if extra:
        msg = (f"{path}: {size} bytes is not a multiple of {RECORD_SIZE}; "
               f"{extra} trailing byte(s)")
        if strict:
            raise ValueError(msg)
        warnings.warn(msg + " ignored", stacklevel=2)
    raw = np.fromfile(path, dtype=RECORD_DTYPE, count=size // RECORD_SIZE)
    if strict and len(raw) != SAMPLES_PER_DWELL:
        raise ValueError(f"{path}: {len(raw)} records, expected {SAMPLES_PER_DWELL}")
    return raw


def read_file(path, dwell=-1, strict=False):
    """Read and decode one ``.dat``.  Returns an ``OUTPUT_DTYPE`` array."""
    return decode(read_raw(path, strict=strict), dwell=dwell)


def find_files(root, pattern="MWB_*.dat"):
    """All burst files under ``root``, ordered by filename (i.e. by time).

    Accepts a single file, a directory, or a glob pattern.
    """
    if os.path.isfile(root):
        return [root]
    if any(ch in root for ch in "*?["):
        found = glob.glob(root, recursive=True)
    else:
        found = glob.glob(os.path.join(root, "**", pattern), recursive=True)
    return sorted(found, key=os.path.basename)


def read_tree(root, pattern="MWB_*.dat", dedupe=True, sort=True, strict=False,
              return_files=False):
    """Read every burst under ``root`` into one flat time series.

    Parameters
    ----------
    root : str
        Directory to search recursively, a single ``.dat``, or a glob.
    dedupe : bool
        Drop records that repeat an already-seen timestamp, keeping the
        first.  Bursts do not overlap in the ARCTERX-2023 data, so this
        normally removes nothing.
    sort : bool
        Sort by time.  Filename order is already chronological, so this is
        a safety net rather than a reordering.
    strict : bool
        Raise instead of warning on short or misaligned files.
    return_files : bool
        Also return the list of files actually read, in ``dwell`` order.

    Returns
    -------
    structured ndarray of ``OUTPUT_DTYPE``, flat and chronological.
    """
    files = find_files(root, pattern)
    if not files:
        raise FileNotFoundError(f"no {pattern} found under {root!r}")

    blocks, used = [], []
    for path in files:
        raw = read_raw(path, strict=strict)
        if len(raw) == 0:
            continue
        blocks.append(decode(raw, dwell=len(used)))
        used.append(path)
    if not blocks:
        raise ValueError(f"no non-empty files under {root!r}")

    ts = np.concatenate(blocks)
    if sort:
        ts = ts[np.argsort(ts["time"], kind="stable")]
    if dedupe:
        keep = np.ones(len(ts), dtype=bool)
        keep[1:] = ts["time"][1:] != ts["time"][:-1]
        ts = ts[keep]
    return (ts, used) if return_files else ts


# ---------------------------------------------------------------- output ---

def save(ts, path):
    """Write a decoded series to ``.npz`` or ``.csv`` (chosen by extension)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".npz":
        np.savez_compressed(path, **{n: ts[n] for n in ts.dtype.names})
    elif ext == ".csv":
        names = [n for n in ts.dtype.names if n != "datetime"]
        with open(path, "w", newline="") as fh:
            fh.write("datetime," + ",".join(names) + "\n")
            for r in ts:
                fh.write(str(r["datetime"]) + "," +
                         ",".join(repr(r[n].item()) for n in names) + "\n")
    else:
        raise ValueError(f"unsupported output extension {ext!r}; use .npz or .csv")


def _summary(ts, files):
    dt = np.diff(ts["time"])
    lines = [
        f"files read            : {len(files)}",
        f"records               : {len(ts)}",
        f"time span             : {ts['datetime'][0]}  ->  {ts['datetime'][-1]}",
        f"sample spacing (ms)   : median {np.median(dt):.0f}, "
        f"min {dt.min()}, max {dt.max()}",
        f"gaps > 1 s            : {int((dt > 1000).sum())}",
        f"lat                   : {ts['lat'].min():.6f} .. {ts['lat'].max():.6f}",
        f"lon                   : {ts['lon'].min():.6f} .. {ts['lon'].max():.6f}",
        f"speed over ground m/s : {ts['sog'].min():.3f} .. {ts['sog'].max():.3f}",
        f"satellites            : {ts['numsats'].min()} .. {ts['numsats'].max()}",
        f"pdop                  : {ts['pdop'].min():.2f} .. {ts['pdop'].max():.2f}",
        f"fix types present     : {sorted(np.unique(ts['fix_type']).tolist())}",
    ]
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("root", help="directory to search, a single .dat, or a glob")
    p.add_argument("-o", "--output", help="write the series to a .npz or .csv")
    p.add_argument("--pattern", default="MWB_*.dat")
    p.add_argument("--strict", action="store_true",
                   help="error out on short or misaligned files")
    p.add_argument("--no-dedupe", action="store_true")
    a = p.parse_args(argv)

    ts, files = read_tree(a.root, pattern=a.pattern, dedupe=not a.no_dedupe,
                          strict=a.strict, return_files=True)
    print(_summary(ts, files))
    if a.output:
        save(ts, a.output)
        print(f"wrote {a.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
