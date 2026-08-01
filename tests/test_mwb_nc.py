#! /usr/bin/env python3

"""Round-trip check: a netCDF written by mwb_nc must reproduce the raw .dat.

Verifies, for each written file:
  1. every variable read back equals the value decoded from the .dat bursts;
  2. the (time, sampling_dwell) layout flattens chronologically the way
     MATLAB's ``val(:)`` does in loadData();
  3. every dwell of the CORDC-delivered file, if one is given, is present.

    python3 test_mwb_nc.py <written.nc> <raw_dir> [delivered.nc]
"""

import os
import sys

import netCDF4
import numpy as np

import mwb_dat
import mwb_nc


def main(new_path, raw_root, delivered=None, start=None, end=None):
    fields, files = mwb_nc.collect_dwells(raw_root, start=start, end=end)
    ds = netCDF4.Dataset(new_path)
    for v in ds.variables.values():
        v.set_auto_mask(False)
        v.set_auto_scale(False)

    ndwell, nsamp = fields["time"].shape
    ok = True

    # ---- 1. dimensions -----------------------------------------------------
    # The burst axis is named sampling_period, not time: CF forbids a
    # multidimensional variable sharing the name of one of its own dimensions.
    DIMS = ("sampling_period", "sampling_dwell")
    assert ds.dimensions[DIMS[0]].size == ndwell, "sampling_period dimension mismatch"
    assert ds.dimensions[DIMS[1]].size == nsamp, "sampling_dwell mismatch"
    assert "time" not in ds.dimensions, \
        "a dimension named 'time' would shadow the 2-D time variable"
    for name in ("time", "u", "v", "lat", "lon"):
        assert ds.variables[name].dimensions == DIMS, \
            f"{name} has dimensions {ds.variables[name].dimensions}"

    # ---- 2. values ---------------------------------------------------------
    checks = [("time", "time"), ("u", "u"), ("v", "v"), ("w", "w"),
              ("lat", "lat"), ("lon", "lon"), ("sog", "sog"), ("cog", "cog"),
              ("numsats", "numsats"), ("pdop", "pdop"), ("height", "height"),
              ("hmsl", "hmsl"), ("h_acc", "h_acc"), ("v_acc", "v_acc"),
              ("s_acc", "s_acc"), ("head_acc", "head_acc"), ("nano", "nano"),
              ("itow", "itow"), ("fix_type", "fix_type")]
    print(f"{'variable':10s} {'max |nc - dat|':>16s}")
    for ncname, myname in checks:
        if ncname not in ds.variables:
            continue
        a = np.asarray(ds.variables[ncname][:], dtype=np.float64)
        b = np.asarray(fields[myname], dtype=np.float64)
        err = np.abs(a - b).max()
        print(f"  {ncname:8s} {err:16.6g}")
        if err != 0.0:
            ok = False
            print(f"    ^^ MISMATCH in {ncname}")

    # ---- 3. MATLAB flattening order ---------------------------------------
    # ncread returns the transpose of the CDL shape, and val(:) is column-major,
    # so MATLAB's flattened vector equals C-order ravel of (time, sampling_dwell).
    flat = np.asarray(ds.variables["time"][:]).ravel(order="C")
    d = np.diff(flat)
    print(f"\nflattened order: {len(flat)} samples, "
          f"monotonic={bool((d > 0).all())}, "
          f"median dt={np.median(d):.0f} ms, max dt={d.max() / 1000:.0f} s")
    if not (d > 0).all():
        ok = False
        print("    ^^ flattened time is not strictly increasing")

    ts = mwb_dat.read_tree(raw_root)
    sel = ts[(ts["time"] >= flat.min()) & (ts["time"] <= flat.max())]
    if len(sel) != len(flat) or not np.array_equal(sel["time"], flat):
        ok = False
        print("    ^^ flattened series differs from mwb_dat.read_tree()")
    else:
        print("flattened series matches mwb_dat.read_tree() exactly")

    # ---- 4. delivered-file coverage ---------------------------------------
    if delivered and os.path.exists(delivered):
        od = netCDF4.Dataset(delivered)
        od.variables["time"].set_auto_mask(False)
        od.variables["time"].set_auto_scale(False)
        ot = od.variables["time"][:].astype(np.int64)
        starts = set(fields["time"][:, 0].tolist())
        missing = [i for i in range(ot.shape[0]) if int(ot[i, 0]) not in starts]
        print(f"\ndelivered file: {ot.shape[0]} dwells, "
              f"{len(missing)} absent from the new file")
        extra = ndwell - (ot.shape[0] - len(missing))
        print(f"new file adds {extra} dwells "
              f"(+{100 * extra / ot.shape[0]:.1f}%)")
        if missing:
            ok = False
            print(f"    ^^ delivered dwells missing: {missing[:10]}")

    print("\n" + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2],
                  sys.argv[3] if len(sys.argv) > 3 else None))
