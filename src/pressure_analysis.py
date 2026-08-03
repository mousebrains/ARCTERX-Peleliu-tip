#! /usr/bin/env python3

"""Regenerate the headline numbers in PRESSURE_ANALYSIS.md.

    python3 src/pressure_analysis.py

`pressure_array.py` was a library with no driver, so every number in
PRESSURE_ANALYSIS.md was asserted rather than computed -- the same failure that
put an unreproducible -0.59 and 5.2x into the drifter half.  This closes that
gap for the tidal results.

It does NOT yet cover the residual-spectrum, noise-floor or gradient-inversion
sections of PRESSURE_ANALYSIS.md; those remain unreproduced.  See the coverage
note at the bottom of the output.
"""

from __future__ import annotations

import argparse
import glob
import os

import numpy as np

import pressure_array as pa
from paths import DATA

# The M2 phase gradient is the one number in PRESSURE_ANALYSIS.md that can be
# checked against another of its own numbers without any fitting at all, so the
# consistency test is built in rather than left to a reader.
DOC_PHASE_SPREAD_DEG = 2.3        # section 4
DOC_M2_GRADIENT_DEG_KM = 0.761    # section 8.1


def sites():
    return [os.path.basename(f).replace("_1min.npz", "")
            for f in sorted(glob.glob(os.path.join(DATA, "pressure", "1min",
                                                   "*_1min.npz")))]


def fit_all(names):
    """Harmonic fit at every gauge.  Returns records and a local x/y frame."""
    out = []
    for s in names:
        t, dep, temp, lat, lon, rate = pa.load(s)
        h = pa.harmonic(t, dep)
        out.append(dict(site=s, lat=lat, lon=lon, amp=h["amp"],
                        var=h["var_explained"], n=int(np.isfinite(dep).sum()),
                        resid=h["residual"]))
    lat = np.array([r["lat"] for r in out])
    lon = np.array([r["lon"] for r in out])
    lat0 = lat.mean()
    x = (lon - lon.mean()) * 111.32 * np.cos(np.radians(lat0))
    y = (lat - lat0) * 110.57
    return out, x, y


def cotidal(recs, x, y, con, period_h):
    """Plane fit to the COMPLEX tidal constant, as PRESSURE_ANALYSIS.md 8.1 says.

    The phase gradient at the array centre is Im(grad Z / Z0); fitting the phase
    angle as a scalar gives the same answer here but is not equivalent in
    general, so the complex form is used.
    """
    a = np.array([r["amp"][con][0] for r in recs])
    p = np.radians(np.array([r["amp"][con][1] for r in recs]))
    Z = a * np.exp(1j * p)
    A = np.column_stack([np.ones_like(x), x, y])
    coef, *_ = np.linalg.lstsq(A, Z, rcond=None)
    resid = Z - A @ coef
    dof = len(x) - 3
    s2 = np.sum(np.abs(resid) ** 2) / dof
    cov = s2 * np.linalg.inv(A.T @ A).real
    Z0, Zx, Zy = coef
    gx = np.degrees((Zx / Z0).imag)
    gy = np.degrees((Zy / Z0).imag)
    g = float(np.hypot(gx, gy))
    sx = np.degrees(np.sqrt(cov[1, 1]) / abs(Z0))
    sy = np.degrees(np.sqrt(cov[2, 2]) / abs(Z0))
    sg = float(np.hypot(sx * abs(gx), sy * abs(gy)) / max(g, 1e-12))
    brg = float(np.degrees(np.arctan2(gx, gy)) % 360)
    speed = 360.0 / (g * period_h * 3600) * 1000 if g > 0 else np.nan
    amp_grad = float(np.hypot((Zx / Z0).real, (Zy / Z0).real) * abs(Z0) * 1000)
    return dict(grad=g, sigma=sg, bearing=brg, speed=speed, amp_grad=amp_grad)


def projected_extent(x, y):
    """Array extent projected onto every direction -- min and max, in km."""
    e = [((x * np.cos(t) + y * np.sin(t)).max()
          - (x * np.cos(t) + y * np.sin(t)).min())
         for t in np.radians(np.arange(0, 180, 1.0))]
    return float(np.min(e)), float(np.max(e))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args()

    names = sites()
    recs, x, y = fit_all(names)

    print(f"PRESSURE ARRAY -- {len(names)} gauges\n")
    print(f"{'site':<6}{'M2 amp m':>10}{'M2 phase':>10}{'var expl':>10}{'n':>9}")
    for r in recs:
        a, p = r["amp"]["M2"]
        print(f"{r['site']:<6}{a:>10.4f}{p:>10.2f}{r['var']:>10.5f}{r['n']:>9d}")

    A = np.array([r["amp"]["M2"][0] for r in recs])
    P = np.array([r["amp"]["M2"][1] for r in recs])
    V = np.array([r["var"] for r in recs])
    spread = float(P.max() - P.min())
    print(f"\nM2 amplitude : mean {A.mean():.4f} m, sd {A.std(ddof=1):.4f}, "
          f"range {A.min():.4f}-{A.max():.4f}")
    print(f"M2 phase     : spread {spread:.2f} deg, sd {P.std(ddof=1):.2f}")
    print(f"variance expl: mean {100*V.mean():.2f} %, min {100*V.min():.2f} %")

    print("\nCO-TIDAL CHART (plane fit to the complex tidal constant)")
    print(f"{'con':<5}{'grad deg/km':>13}{'sigma':>9}{'sig':>7}"
          f"{'bearing':>9}{'speed m/s':>11}")
    per = {"M2": 12.4206012, "S2": 12.0, "O1": 25.8193417, "K1": 23.9344697}
    res = {}
    for c, ph in per.items():
        r = cotidal(recs, x, y, c, ph)
        res[c] = r
        sig = r["grad"] / r["sigma"] if r["sigma"] > 0 else 0.0
        print(f"{c:<5}{r['grad']:>13.3f}{r['sigma']:>9.3f}{sig:>7.1f}"
              f"{r['bearing']:>9.0f}{r['speed']:>11.1f}")
    print(f"M2 amplitude gradient: {res['M2']['amp_grad']:.2f} mm/km")

    lo, hi = projected_extent(x, y)
    g = res["M2"]["grad"]
    print(f"\nCONSISTENCY: array projects onto {lo:.1f}-{hi:.1f} km.")
    print(f"  measured gradient {g:.3f} deg/km -> phase spread "
          f"{g*lo:.1f}-{g*hi:.1f} deg; measured spread is {spread:.2f} deg.")
    print(f"  PRESSURE_ANALYSIS.md 8.1 quotes {DOC_M2_GRADIENT_DEG_KM} deg/km, "
          f"which requires a spread of "
          f"{DOC_M2_GRADIENT_DEG_KM*lo:.1f}-{DOC_M2_GRADIENT_DEG_KM*hi:.1f} deg")
    print(f"  against the {DOC_PHASE_SPREAD_DEG} deg its own section 4 reports.")
    ok = DOC_M2_GRADIENT_DEG_KM * lo <= DOC_PHASE_SPREAD_DEG
    print(f"  -> documented gradient and documented spread are "
          f"{'CONSISTENT' if ok else 'INCONSISTENT'}.")

    c_meas = 360.0 / (g * 12.4206012 * 3600) * 1000
    c_doc = 360.0 / (DOC_M2_GRADIENT_DEG_KM * 12.4206012 * 3600) * 1000
    print(f"\n  phase speed: measured {c_meas:.1f} m/s (implied depth "
          f"{c_meas**2/9.81:.0f} m); documented {c_doc:.1f} m/s "
          f"(implied depth {c_doc**2/9.81:.0f} m).")

    print("\nNOT COVERED HERE (still unreproduced in PRESSURE_ANALYSIS.md):")
    print("  residual spectra and the 2.2 cm / 3.3 cm noise-floor split (6),")
    print("  the 83 % non-barotropic fraction and 2-7x current bias (4),")
    print("  the C05 rotation and subtidal regime classification (4, 8).")


if __name__ == "__main__":
    main()
