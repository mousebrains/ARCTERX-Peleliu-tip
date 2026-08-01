#! /usr/bin/env python3

"""Write CORDC-format GPS timeseries netCDF files from raw mwb ``.dat`` bursts.

Produces one ``mwb<serial>d<NN>_gps_timeseries.nc`` per drifter, spanning that
drifter's whole time range, in a CF-1.13 / ACDD-1.3 conformant version of the
CORDC-delivered layout: every variable is ``(sampling_period, sampling_dwell)``, i.e.
blocks of 2048 samples per wave-processing period.  MATLAB's ``ncread`` returns
those transposed as (2048, ndwell), so ``val(:)`` flattens to a chronological
column exactly as ``loadData()`` in ``explore_drifter_paths.m`` expects.

No trimming or quality filtering is applied -- every burst found is written.

DIFFERENCES FROM THE CORDC-DELIVERED FILES
------------------------------------------
By default this writer keeps the raw precision that the delivered files threw
away, and adds the GPS fields they omit entirely:

  * ``pdop`` is float32 with 0.01 resolution.  The delivered files store
    ``round(pDOP/100)`` in a signed byte with ``scale_factor = 1``, so their
    pdop only ever takes the values 1 or 2 while the true range is 0.30-2.52.
  * ``sog`` / ``cog`` are float32.  The delivered files quantize to 0.1 m/s
    and 0.1 deg via ``short`` with ``scale_factor = 0.1``.
  * Added: ``height``, ``hmsl``, ``h_acc``, ``v_acc``, ``s_acc``, ``head_acc``,
    ``t_acc``, ``nano``, ``itow``, ``fix_type``, ``flags``, ``flags2``,
    ``valid``, and ``source_file`` (the ``.dat`` each dwell came from).

``loadData()`` reads a fixed list of ten variable names, so the extra
variables are inert there, and ``ncread`` applies ``scale_factor``
automatically -- both changes are transparent to the existing pipeline.

Pass ``--cordc-compat`` to reproduce the delivered encodings exactly instead
(short sog/cog, byte pdop, core variables only).

USAGE
-----
    python3 mwb_nc.py --input  /Volumes/SeaChest/ARCTERX/2023/Wake/MiniWaveBuoys \\
                      --output ./netcdf.new \\
                      --serials 458,788,790,793 \\
                      --start 2023-05-22 --end 2023-05-23 \\
                      --deployment 458=2,788=1,790=1,793=2

``--input`` may be a directory of per-serial subdirectories (``458/``,
``788/``, ...), or a single drifter's directory.  Bursts are found
recursively as ``<dir>/**/MWB_*.dat``.

Requires numpy, plus either ``netCDF4`` (preferred: writes netCDF-4/HDF5 with
compression, matching the delivered files) or ``scipy`` (falls back to
netCDF-3 classic, which MATLAB reads equally well).
"""

from __future__ import annotations

import argparse
import datetime as _dt
import os
import re
import sys
import warnings

import numpy as np

import mwb_dat

__all__ = ["write_netcdf", "collect_dwells", "main"]

FILL_F4 = np.float32(9.96921e36)
FILL_F8 = -2147483647.0
FILL_I2 = np.int16(-32767)
FILL_I4 = np.int32(-2147483647)
FILL_I1 = np.int8(-127)

_DEPLOY_STRLEN = 56
_FNAME_STRLEN = 64

#: Attached to every geophysical variable so CF can locate it in space, time
#: and trajectory, and resolve the horizontal datum.
COORDS = {"coordinates": "time lon lat deployment", "grid_mapping": "wgs84"}

#: UBX-NAV-PVT fixType enumeration (u-blox interface description).  Only value
#: 3 occurs in the ARCTERX-2023 record.
FIX_TYPE_VALUES = np.array([0, 1, 2, 3, 4, 5], dtype=np.int8)
FIX_TYPE_MEANINGS = ("no_fix dead_reckoning_only 2D_fix 3D_fix "
                     "GNSS_plus_dead_reckoning time_only_fix")


# ------------------------------------------------------- backend adapter ---

class _Writer:
    """Thin adapter over netCDF4 or scipy.io.netcdf_file."""

    def __init__(self, path, compress=True, complevel=4):
        self.path = path
        self._complevel = complevel
        try:
            import netCDF4
            self._kind = "netcdf4"
            self._ds = netCDF4.Dataset(path, "w", format="NETCDF4")
            self._compress = compress
        except ImportError:
            try:
                from scipy.io import netcdf_file
            except ImportError:
                raise ImportError(
                    "writing netCDF needs either the 'netCDF4' or the 'scipy' "
                    "package; install one with e.g. `pip install netCDF4`") from None
            self._kind = "scipy"
            self._ds = netcdf_file(path, "w")
            self._compress = False

    @property
    def kind(self):
        return self._kind

    def dim(self, name, size):
        self._ds.createDimension(name, size)

    def var(self, name, dtype, dims, attrs=None, data=None, fill=None):
        attrs = dict(attrs or {})
        if self._kind == "netcdf4":
            kw = {}
            if fill is not None:
                kw["fill_value"] = fill
            if self._compress and dims:
                kw.update(zlib=True, complevel=self._complevel)
            v = self._ds.createVariable(name, dtype, dims, **kw)
            v.set_auto_maskandscale(False)   # write pre-packed values verbatim
            for k, val in attrs.items():
                v.setncattr(k, val)
        else:
            code = np.dtype(dtype).char
            code = {"b": "b", "B": "b", "h": "h", "i": "i", "l": "i",
                    "f": "f", "d": "d", "S": "c"}.get(code, code)
            v = self._ds.createVariable(name, code, dims)
            if fill is not None:
                v._FillValue = fill
            for k, val in attrs.items():
                setattr(v, k, val)
        if data is not None:
            with warnings.catch_warnings():
                # netCDF4 <= 1.7.4 sets .shape on an ndarray internally, which
                # numpy >= 2.5 deprecates.  Library-internal, harmless here.
                warnings.filterwarnings("ignore", category=DeprecationWarning,
                                        message=".*shape on a NumPy array.*")
                v[:] = data
        return v

    def globals(self, attrs):
        for k, val in attrs.items():
            if self._kind == "netcdf4":
                self._ds.setncattr(k, val)
            else:
                setattr(self._ds, k, val)

    def close(self):
        self._ds.close()


def _iso_duration(seconds):
    """Seconds -> an ISO 8601 duration, as ACDD wants for time_coverage_*."""
    if not np.isfinite(seconds):
        return "unknown"
    days, rem = divmod(float(seconds), 86400.0)
    hours, rem = divmod(rem, 3600.0)
    mins, secs = divmod(rem, 60.0)
    out = "P" + (f"{int(days)}D" if days else "")
    tail = ((f"{int(hours)}H" if hours else "")
            + (f"{int(mins)}M" if mins else "")
            + (f"{secs:g}S" if secs else ""))
    return out + ("T" + tail if tail else "") if (days or tail) else "PT0S"


def _chararray(strings, strlen):
    """List of str -> (n, strlen) array of single bytes, NUL padded."""
    out = np.zeros((len(strings), strlen), dtype="S1")
    for i, s in enumerate(strings):
        b = s.encode("utf-8")[:strlen]
        out[i, :len(b)] = np.frombuffer(b, dtype="S1")
    return out


# --------------------------------------------------------- data assembly ---

def _stamp_ms(path):
    """UTC ms from the ``MWB_<yyyymmddHHMMSS>.dat`` burst-command stamp."""
    m = re.search(r"MWB_(\d{14})", os.path.basename(path))
    if not m:
        return None
    d = _dt.datetime.strptime(m.group(1), "%Y%m%d%H%M%S")
    return int(d.replace(tzinfo=_dt.timezone.utc).timestamp() * 1000)


def collect_dwells(root, pattern="MWB_*.dat", start=None, end=None):
    """Read every burst under ``root`` into (ndwell, 2048) arrays.

    A burst is kept if any of its records falls within [start, end]; bursts are
    never split, so the 2048-sample block structure is preserved.  Bursts are
    ordered by time.  Short bursts (none in the ARCTERX-2023 data) are padded
    with fill values.

    Returns ``(fields, files)`` where ``fields`` maps ``OUTPUT_DTYPE`` names to
    (ndwell, nsamp) arrays and ``files`` is the source path of each dwell.
    """
    paths = mwb_dat.find_files(root, pattern)
    if not paths:
        raise FileNotFoundError(f"no {pattern} found under {root!r}")
    if start is not None or end is not None:
        # Cheap prefilter on the burst-command stamp in the filename, so we
        # never open files from unrelated dates.  The stamp leads the first
        # record by a few seconds and a burst lasts at most ~35 min, so an
        # hour of margin either side cannot drop a relevant file.
        paths = [p for p in paths
                 if _stamp_ms(p) is None
                 or ((start is None or _stamp_ms(p) >= start - 3_600_000)
                     and (end is None or _stamp_ms(p) <= end + 3_600_000))]
        if not paths:
            raise ValueError(f"no bursts under {root!r} within the requested range")

    blocks, files = [], []
    for p in paths:
        raw = mwb_dat.read_raw(p)
        if len(raw) == 0:
            continue
        rec = mwb_dat.decode(raw, dwell=len(files))
        if start is not None and rec["time"].max() < start:
            continue
        if end is not None and rec["time"].min() > end:
            continue
        blocks.append(rec)
        files.append(p)
    if not blocks:
        raise ValueError(f"no bursts under {root!r} within the requested range")

    order = np.argsort([b["time"][0] for b in blocks], kind="stable")
    blocks = [blocks[i] for i in order]
    files = [files[i] for i in order]

    nsamp = max(len(b) for b in blocks)
    if nsamp != mwb_dat.SAMPLES_PER_DWELL:
        warnings.warn(f"longest burst is {nsamp} samples, expected "
                      f"{mwb_dat.SAMPLES_PER_DWELL}", stacklevel=2)

    fields = {}
    for name in mwb_dat.OUTPUT_DTYPE.names:
        if name in ("datetime", "dwell"):
            continue
        col = np.empty((len(blocks), nsamp), dtype=mwb_dat.OUTPUT_DTYPE[name])
        col[:] = 0
        for i, b in enumerate(blocks):
            col[i, :len(b)] = b[name]
        fields[name] = col
    fields["_valid"] = np.array(
        [[j < len(b) for j in range(nsamp)] for b in blocks], dtype=bool)
    return fields, files


def _sampling_rate(time_ms, valid):
    """Nominal Hz per dwell, from the median sample spacing."""
    out = np.zeros(time_ms.shape[0], dtype=np.int8)
    for i in range(time_ms.shape[0]):
        t = time_ms[i][valid[i]]
        out[i] = 0 if len(t) < 2 else int(round(1000.0 / np.median(np.diff(t))))
    return out


# ---------------------------------------------------------------- writer ---

def write_netcdf(path, fields, files, *, serial, deployment,
                 deployment_name="ARCTERX May 2023", cordc_compat=False,
                 title=None, summary=None, compress=True, complevel=4,
                 publisher=None):
    """Write one drifter's bursts to a CORDC-format netCDF file."""
    ndwell, nsamp = fields["time"].shape
    valid = fields["_valid"]
    buoy_id = f"mwb{serial}d{deployment:02d}"

    t = fields["time"].astype(np.float64)
    t[~valid] = FILL_F8
    rate = _sampling_rate(fields["time"], valid)

    def f4(name, fill=FILL_F4):
        a = fields[name].astype(np.float32).copy()
        a[~valid] = fill
        return a

    lat, lon = f4("lat"), f4("lon")
    fin = valid & np.isfinite(fields["lat"]) & np.isfinite(fields["lon"])

    w = _Writer(path, compress=compress, complevel=complevel)
    # CF forbids a multidimensional variable sharing the name of one of its own
    # dimensions (it would be mistaken for a coordinate variable), so the burst
    # axis is named sampling_period rather than time.  Dimension names are not
    # referenced by ncread(), so this is invisible to the MATLAB pipeline.
    w.dim("sampling_period", ndwell)
    w.dim("sampling_dwell", nsamp)
    w.dim("deploy_strlen", _DEPLOY_STRLEN)
    DIMS = ("sampling_period", "sampling_dwell")

    w.var("deployment", "S1", ("deploy_strlen",),
          {"standard_name": "platform_name", "long_name": "deployment",
           "cf_role": "trajectory_id", "coverage_content_type": "referenceInformation"},
          data=_chararray([deployment_name], _DEPLOY_STRLEN)[0])
    w.var("sampling_dwell", "i2", (),
          {"long_name": "number of samples collected per wave processing period",
           "units": "1", "coverage_content_type": "auxiliaryInformation"},
          data=np.int16(nsamp), fill=np.int16(0))
    w.var("sampling_rate", "i1", ("sampling_period",),
          {"units": "Hz", "long_name": "GPS sampling frequency",
           "coverage_content_type": "auxiliaryInformation",
           "comment": "Defines the sampling frequency for GPS velocity measurements"},
          data=rate, fill=np.int8(0))
    # calendar "standard", not "utc": these values are a linear count of
    # milliseconds that ignores leap seconds (CF-1.13 section 4.4.1), which is
    # what the receiver's UTC calendar fields encode and what POSIX readers,
    # including MATLAB's datetime(...,ConvertFrom="posixtime"), expect.
    # "gregorian" is a deprecated alias for "standard" as of CF-1.9.
    w.var("time", "f8", DIMS,
          {"standard_name": "time", "long_name": "time of GPS measurement",
           "units": "milliseconds since 1970-01-01", "calendar": "standard",
           "axis": "T", "coverage_content_type": "coordinate"},
          data=t, fill=FILL_F8)

    for name, sn, ln in (("u", "eastward_sea_water_velocity", "GPS velE measurement"),
                         ("v", "northward_sea_water_velocity", "GPS velN measurement"),
                         ("w", "upward_sea_water_velocity", "GPS velZ measurement")):
        w.var(name, "f4", DIMS,
              {"standard_name": sn, "units": "m s-1", "long_name": ln,
               "coverage_content_type": "physicalMeasurement", **COORDS},
              data=f4(name, np.float32(-2.147484e9)), fill=np.float32(-2.147484e9))

    w.var("lat", "f4", DIMS,
          {"standard_name": "latitude", "long_name": "latitude",
           "units": "degrees_north", "axis": "Y",
           "coverage_content_type": "coordinate"},
          data=lat, fill=FILL_F4)
    w.var("lon", "f4", DIMS,
          {"standard_name": "longitude", "long_name": "longitude",
           "units": "degrees_east", "axis": "X",
           "coverage_content_type": "coordinate"},
          data=lon, fill=FILL_F4)
    w.var("wgs84", "i1", (),
          {"grid_mapping_name": "latitude_longitude",
           "longitude_of_prime_meridian": np.float32(0.0),
           "semi_major_axis": np.float32(6378137.0),
           "inverse_flattening": 298.257223563},
          data=np.int8(-127), fill=np.int8(-127))

    coords = dict(COORDS)
    if cordc_compat:
        for name, sn, ln, units in (
                ("sog", "platform_speed_wrt_ground",
                 "GPS speed over ground measurement", "m s-1"),
                ("cog", "platform_course",
                 "GPS course over ground measurement", "degree")):
            packed = np.floor(fields[name].astype(np.float64) * 10 + 0.5)
            packed = np.where(valid, packed, FILL_I2).astype(np.int16)
            w.var(name, "i2", DIMS,
                  {"scale_factor": np.float32(0.1), "standard_name": sn,
                   "long_name": ln, "units": units,
                   "coverage_content_type": "physicalMeasurement", **coords},
                  data=packed, fill=FILL_I2)
        pd = np.where(valid, np.round(fields["pdop"]), FILL_I1).astype(np.int8)
        w.var("pdop", "i1", DIMS,
              {"scale_factor": 1.0, "long_name": "GPS dilution of precision",
               "units": "dB", "coverage_content_type": "qualityInformation",
               "comment": "quality indicator for GPS measurements", **coords},
              data=pd, fill=FILL_I1)
    else:
        w.var("sog", "f4", DIMS,
              {"standard_name": "platform_speed_wrt_ground",
               "long_name": "GPS speed over ground measurement",
               "units": "m s-1",
               "coverage_content_type": "physicalMeasurement", **coords},
              data=f4("sog"), fill=FILL_F4)
        w.var("cog", "f4", DIMS,
              {"standard_name": "platform_course",
               "long_name": "GPS course over ground measurement",
               "units": "degree",
               "coverage_content_type": "physicalMeasurement", **coords},
              data=f4("cog"), fill=FILL_F4)
        w.var("pdop", "f4", DIMS,
              {"long_name": "GPS position dilution of precision", "units": "1",
               "coverage_content_type": "qualityInformation",
               "comment": "quality indicator for GPS measurements; "
                          "0.01 resolution as recorded by the receiver", **coords},
              data=f4("pdop"), fill=FILL_F4)

    ns = np.where(valid, fields["numsats"], 127).astype(np.int8)
    w.var("numsats", "i1", DIMS,
          {"long_name": "GPS satellites used in the navigation solution",
           "units": "1", "coverage_content_type": "qualityInformation",
           "comment": "quality indicator for GPS measurements", **coords},
          data=ns, fill=np.int8(127))

    if not cordc_compat:
        # standard_name is set only where a real CF standard name exists.
        # GPS accuracy estimates, DOP, satellite counts and receiver bitfields
        # have none in the CF table (checked against v94), and inventing one
        # would break CF compliance to satisfy an ACDD recommendation.
        # CF 4.3: a variable carrying a vertical-coordinate standard name must
        # also declare `positive`.  Both of these are heights measured upward
        # from their respective datums, so both get positive = "up".  Omitting
        # it on hmsl (standard_name "altitude") is a hard CF error, caught by
        # the IOOS compliance-checker at cf:1.11 but not by cfchecker at 1.8.
        for name, units, ln, sn, cct, pos in (
                ("height", "m", "height above the WGS-84 ellipsoid",
                 "height_above_reference_ellipsoid", "physicalMeasurement", "up"),
                ("hmsl", "m", "height above mean sea level",
                 "altitude", "physicalMeasurement", "up"),
                ("h_acc", "m", "GPS horizontal accuracy estimate",
                 None, "qualityInformation", None),
                ("v_acc", "m", "GPS vertical accuracy estimate",
                 None, "qualityInformation", None),
                ("s_acc", "m s-1", "GPS speed accuracy estimate",
                 None, "qualityInformation", None),
                ("head_acc", "degree", "GPS heading accuracy estimate",
                 None, "qualityInformation", None)):
            attrs = {"units": units, "long_name": ln,
                     "coverage_content_type": cct, **coords}
            if sn:
                attrs["standard_name"] = sn
            if pos:
                attrs["positive"] = pos
            w.var(name, "f4", DIMS, attrs, data=f4(name), fill=FILL_F4)
        for name, dtype, fill, units, ln in (
                ("t_acc", "i4", FILL_I4, "ns", "GPS time accuracy estimate"),
                ("nano", "i4", FILL_I4, "ns",
                 "sub-second residual of the UTC time stamp"),
                ("itow", "i4", FILL_I4, "ms", "GPS time of week")):
            a = np.where(valid, fields[name], fill).astype(dtype)
            w.var(name, dtype, DIMS,
                  {"units": units, "long_name": ln,
                   "coverage_content_type": "auxiliaryInformation", **coords},
                  data=a, fill=fill)
        for name, ln, cmt in (
                ("fix_type", "GNSS fix type", None),
                ("flags", "UBX-NAV-PVT fix status bitfield",
                 "bit 0 gnssFixOK, bit 1 diffSoln, bits 2-4 psmState, "
                 "bit 5 headVehValid, bits 6-7 carrSoln"),
                ("flags2", "UBX-NAV-PVT additional flags bitfield",
                 "bit 5 confirmedAvai, bit 6 confirmedDate, bit 7 confirmedTime"),
                ("valid", "UBX-NAV-PVT time validity bitfield",
                 "bit 0 validDate, bit 1 validTime, bit 2 fullyResolved, "
                 "bit 3 validMag")):
            a = np.where(valid, fields[name], FILL_I1).astype(np.int8)
            attrs = {"long_name": ln, "units": "1",
                     "coverage_content_type": "qualityInformation", **coords}
            if name == "fix_type":
                # Enumerated, so flag_values applies.  The other three mix
                # single- and multi-bit subfields, which flag_masks cannot
                # express faithfully, so they carry a comment instead.
                attrs["flag_values"] = FIX_TYPE_VALUES
                attrs["flag_meanings"] = FIX_TYPE_MEANINGS
            if cmt:
                attrs["comment"] = cmt
            w.var(name, "i1", DIMS, attrs, data=a, fill=FILL_I1)
        w.dim("filename_strlen", _FNAME_STRLEN)
        w.var("source_file", "S1", ("sampling_period", "filename_strlen"),
              {"long_name": "raw .dat burst file each dwell was decoded from",
               "coverage_content_type": "auxiliaryInformation"},
              data=_chararray([os.path.basename(f) for f in files], _FNAME_STRLEN))

    tmin = np.datetime64(int(fields["time"][valid].min()), "ms")
    tmax = np.datetime64(int(fields["time"][valid].max()), "ms")
    now = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    rates = sorted({int(r) for r in rate if r})
    rate_txt = "/".join(str(r) for r in rates) or "unknown"

    latmin, latmax = float(fields["lat"][fin].min()), float(fields["lat"][fin].max())
    lonmin, lonmax = float(fields["lon"][fin].min()), float(fields["lon"][fin].max())
    bounds_wkt = ("POLYGON (({0:.6f} {2:.6f}, {1:.6f} {2:.6f}, {1:.6f} {3:.6f}, "
                  "{0:.6f} {3:.6f}, {0:.6f} {2:.6f}))"
                  ).format(lonmin, lonmax, latmin, latmax)
    dts = np.diff(fields["time"], axis=1)[valid[:, 1:]]
    res_s = float(np.median(dts)) / 1000.0 if dts.size else float("nan")

    w.globals({
        "id": buoy_id,
        "title": ("Near real-time in situ directional wave and sea surface "
                  "temperature\nmeasurements collected in Peleliu Wake by a CORDC "
                  f"Miniature Wave Buoy\n(unit {serial}, deployment {deployment})"),
        "summary": (summary or
                    f"CORDC Miniature Wave Buoy {serial} (deployment {deployment}) was "
                    f"deployed in Peleliu Wake.  In situ GPS velocity and location data\n"
                    f"({ndwell} samples) were decoded from the raw onboard .dat burst "
                    f"files and cover\n{tmin} to {tmax}.\n"
                    f"{nsamp} samples are collected each sampling period; the sampling "
                    f"duration is a function of the sampling frequency.\n"
                    f"This dataset contains GPS records sampled at {rate_txt} Hz.\n"
                    "No trimming or quality-control filtering has been applied: every "
                    "burst present in the raw record is included."),
        "keywords": ("Earth Science, Sea Surface, Oceans, Ocean Waves, Gravity Waves, "
                     "Sea State, Ocean, Pacific Ocean, Western Pacific Ocean, "
                     "Micronesia, Palau, Peleliu Wake"),
        "source": "in situ sea surface measurements",
        "platform": "In Situ Ocean-based Platforms, Buoys, 15 inch hull",
        "instrument": ("Earth Remote Sensing Instruments, Passive Remote Sensing, "
                       "Positioning/Navigation, GPS, GPS Receivers"),
        "institution": ("Coastal Observing Research and Development Center, Scripps "
                        "Institution of Oceanography"),
        "processing_level": ("Raw onboard GPS bursts decoded from the instrument .dat "
                             "files with no trimming, gap filling or quality-control "
                             "filtering applied"),
        "time_coverage_start": str(tmin) + "Z",
        "time_coverage_end": str(tmax) + "Z",
        "time_coverage_duration": _iso_duration(
            (int(fields["time"][valid].max()) - int(fields["time"][valid].min())) / 1000.0),
        "time_coverage_resolution": _iso_duration(res_s),
        "geospatial_lat_min": np.float32(latmin),
        "geospatial_lat_max": np.float32(latmax),
        "geospatial_lat_units": "degrees_north",
        "geospatial_lon_min": np.float32(lonmin),
        "geospatial_lon_max": np.float32(lonmax),
        "geospatial_lon_units": "degrees_east",
        "geospatial_bounds": bounds_wkt,
        "geospatial_bounds_crs": "EPSG:4326",
        "date_created": now,
        "date_modified": now,
        "creator_name": "Coastal Observing Research and Development Center",
        "creator_type": "group",
        "creator_email": "cordc.wavebuoy@sio.ucsd.edu",
        "creator_url": "http://cordc.ucsd.edu/",
        "project": "ARCTERX 2023 Wake experiment, Peleliu (Palau)",
        "license": ("Not specified by the data originator; contact "
                    "cordc.wavebuoy@sio.ucsd.edu for terms of use."),
        "references": "Available upon request",
        "acknowledgement": ("Data collected by the Coastal Observing Research and "
                            "Development Center, Scripps Institution of Oceanography."),
        "comment": ("Decoded directly from the instrument's raw 85-byte UBX-NAV-PVT "
                    "burst records, so this file carries every burst the buoy "
                    "recorded, at the receiver's full precision.  Geospatial vertical "
                    "extent is deliberately not advertised: the GPS height of a "
                    "surface-following buoy is dominated by receiver noise and does "
                    "not describe a sampled depth range."),
        "history": f"{now} written by mwb_nc.py from {len(files)} raw .dat bursts",
        "source": "in situ sea surface measurements",
        "standard_name_vocabulary": "CF Standard Name Table v94",
        "Conventions": "CF-1.13, ACDD-1.3",
        "featureType": "trajectory",
        "cdm_data_type": "trajectory",
        "naming_authority": "edu.ucsd.cordc",
        "keywords_vocabulary": ("Global Change Master Directory (GCMD) Keywords, "
                                "Version 8.5"),
        "deployment": deployment_name,
        # ACDD recommends publisher_*, but these files are derived locally
        # from the raw instrument record; naming a publisher that has not
        # actually published them would be fabricated provenance.  Supply
        # --publisher-name/-email/-url when you do distribute them.
        **{f"publisher_{k}": v for k, v in (publisher or {}).items() if v},
    })
    w.close()
    return {"id": buoy_id, "path": path, "ndwell": ndwell, "nsamp": nsamp,
            "start": tmin, "end": tmax, "backend": w.kind, "nfiles": len(files)}


# ------------------------------------------------------------------- CLI ---

def _parse_day(s, end=False):
    for fmt in ("%Y-%m-%d", "%Y%m%d", "%Y-%m-%dT%H:%M:%S"):
        try:
            d = _dt.datetime.strptime(s, fmt).replace(tzinfo=_dt.timezone.utc)
            if end and fmt in ("%Y-%m-%d", "%Y%m%d"):
                d += _dt.timedelta(days=1) - _dt.timedelta(milliseconds=1)
            return int(d.timestamp() * 1000)
        except ValueError:
            pass
    raise argparse.ArgumentTypeError(f"cannot parse date {s!r}")


def _find_serial_dirs(root, serials=None):
    """Map serial -> directory.  Handles both a parent of serial dirs and a
    single drifter directory."""
    subs = {}
    if os.path.isdir(root):
        for name in sorted(os.listdir(root)):
            p = os.path.join(root, name)
            if os.path.isdir(p) and re.fullmatch(r"\d+", name):
                subs[name] = p
    if not subs:
        m = re.search(r"(\d+)\D*$", os.path.basename(os.path.normpath(root)))
        subs = {m.group(1) if m else "unknown": root}
    if serials:
        missing = [s for s in serials if s not in subs]
        if missing:
            raise SystemExit(f"no directory for serial(s) {', '.join(missing)} "
                             f"under {root}; found {', '.join(subs) or 'none'}")
        subs = {s: subs[s] for s in serials}
    return subs


def _deployment_map(spec, serials, reference_dirs):
    """Resolve serial -> deployment number, from ``spec`` then from any
    existing ``mwb<serial>d<NN>_gps_timeseries.nc`` in ``reference_dirs``."""
    out = {}
    for item in filter(None, (spec or "").split(",")):
        k, _, v = item.partition("=")
        out[k.strip()] = int(v)
    for s in serials:
        if s in out:
            continue
        found = set()
        for d in reference_dirs:
            if not os.path.isdir(d):
                continue
            for fn in os.listdir(d):
                m = re.fullmatch(rf"mwb{s}d(\d+)_gps_timeseries\.nc", fn)
                if m:
                    found.add(int(m.group(1)))
        if len(found) == 1:
            out[s] = found.pop()
        elif len(found) > 1:
            raise SystemExit(
                f"serial {s}: several deployments found ({sorted(found)}); "
                f"disambiguate with --deployment {s}=NN")
        else:
            raise SystemExit(
                f"serial {s}: cannot infer the deployment number; "
                f"pass it with --deployment {s}=NN")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", required=True,
                   help="directory of per-serial subdirectories, or one drifter's dir")
    p.add_argument("-o", "--output", required=True,
                   help="directory to write mwb<serial>d<NN>_gps_timeseries.nc into")
    p.add_argument("--serials", help="comma-separated serials (default: all found)")
    p.add_argument("--deployment", help="deployment numbers, e.g. 458=2,788=1")
    p.add_argument("--start", help="drop bursts entirely before this date (UTC)")
    p.add_argument("--end", help="drop bursts entirely after this date (UTC)")
    p.add_argument("--deployment-name", default="ARCTERX May 2023")
    p.add_argument("--reference", action="append", default=[],
                   help="dir of existing CORDC files, to infer deployment numbers "
                        "(default: --output and the working directory)")
    p.add_argument("--cordc-compat", action="store_true",
                   help="reproduce the delivered encodings exactly (lossy sog/cog/"
                        "pdop, core variables only)")
    p.add_argument("--no-compress", action="store_true")
    p.add_argument("--complevel", type=int, default=4,
                   help="zlib deflate level 1-9 (default: 4)")
    p.add_argument("--publisher-name")
    p.add_argument("--publisher-email")
    p.add_argument("--publisher-url")
    p.add_argument("--pattern", default="MWB_*.dat")
    a = p.parse_args(argv)

    start = _parse_day(a.start) if a.start else None
    end = _parse_day(a.end, end=True) if a.end else None
    serials = [s.strip() for s in a.serials.split(",")] if a.serials else None
    dirs = _find_serial_dirs(a.input, serials)
    refs = a.reference or [a.output, os.getcwd()]
    depl = _deployment_map(a.deployment, list(dirs), refs)

    os.makedirs(a.output, exist_ok=True)
    for serial, src in dirs.items():
        fields, files = collect_dwells(src, a.pattern, start=start, end=end)
        out = os.path.join(a.output,
                           f"mwb{serial}d{depl[serial]:02d}_gps_timeseries.nc")
        info = write_netcdf(out, fields, files, serial=serial,
                            deployment=depl[serial],
                            deployment_name=a.deployment_name,
                            cordc_compat=a.cordc_compat,
                            compress=not a.no_compress, complevel=a.complevel,
                            publisher={'name': a.publisher_name,
                                       'email': a.publisher_email,
                                       'url': a.publisher_url})
        print(f"{info['id']}: {info['ndwell']} dwells x {info['nsamp']} samples "
              f"({info['ndwell'] * info['nsamp']} records) "
              f"{info['start']} -> {info['end']}")
        print(f"  {out}  [{info['backend']}, "
              f"{os.path.getsize(out) / 1e6:.1f} MB]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
