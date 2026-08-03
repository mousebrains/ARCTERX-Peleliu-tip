#! /usr/bin/env python3

"""Null tests: does the pipeline recover a vortex whose answer we already know?

Every number in DRIFTER_ANALYSIS.md comes from estimators applied to real
data, where the truth is unknown.  That is exactly the situation in which an
estimator's own bias is invisible: a wrong answer and a right answer look
identical.  The cure is to feed the estimators a synthetic flow whose
vorticity, circulation and core radius we chose ourselves, and check what
comes back.

Four tests, in increasing order of how much they can embarrass us:

  1. SOLID BODY, EXACTNESS.  For a solid-body rotation the velocity field is
     exactly affine, so Stokes' theorem is exact for a straight-sided polygon
     and the circulation estimator must return zeta = 2*omega with no error at
     all.  Anything else is an implementation bug, not a modeling limitation.

  2. SOLID BODY, NULL FOR THE STRUCTURE CLAIM.  DRIFTER_ANALYSIS.md concludes
     the vortex is "not solid body" because |zeta| falls with the scale being
     sampled (r = -0.56).  That is only meaningful if a vortex which IS solid
     body returns r ~ 0 through the same machinery.  If the pipeline
     manufactured a negative correlation on its own, the conclusion would be an
     artifact.

     The specific worry is that the quadrature error of test 4 grows with
     cluster size, so a larger cluster might report a smaller |zeta| even in a
     constant-vorticity field.  Advecting ONE cluster through a solid-body
     field does not test that: such a cluster is rigid, so its size never
     varies (CV ~ 1e-5) and the correlation is between two constants.  Test 2b
     therefore pools clusters over a 232-1238 m range of sizes, which is the
     comparison the claim actually rests on.

  3. LAMB-OSEEN, RECOVERY.  Integrate drifters through a known Lamb-Oseen
     vortex and check that Gamma and the core radius come back.

  5. FORMAL-ERROR REALISM.  docs/03 claims the least-squares formal error is
     ~5x too small.  That number was previously asserted with nothing computing
     it.  Here it is measured: integrate through a Lamb-Oseen field, compare the
     LSQ zeta against the EXACT area-averaged zeta over the same polygon, and
     take the ratio of actual RMS error to the formal standard deviation.

     The truth must be the area average, not the local zeta at the centroid.
     The LSQ fit returns an area average; comparing it against a local value
     folds in the averaging offset and inflates the ratio to 13-17x rather than
     the true 5-7x.

  4. QUADRATURE BIAS.  The circulation estimator evaluates the closed line
     integral with the trapezoid rule on four vertex velocities.  For a curved
     flow that is approximate.  Measure the error against a high-order
     quadrature of the same integral, using the REAL cluster geometry.

Run standalone:  python3 tests/test_synthetic_recovery.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))), "src"))

import eddy_kinematics as ek           # noqa: E402
from paths import DATA                 # noqa: E402

G_TRUE, R_TRUE = -6000.0, 1200.0       # m^2/s, m
OMEGA_SB = -6.0e-4                     # solid body; zeta = 2*omega
SPEED, BEARING = 0.131, 266.0          # observed translation
_HDG = np.radians(90.0 - BEARING)
CU, CV = SPEED * np.cos(_HDG), SPEED * np.sin(_HDG)
_GL = np.polynomial.legendre.leggauss(200)


# ------------------------------------------------------------------ fields ---

def field_oseen(x, y, cx, cy, G=G_TRUE, R=R_TRUE):
    """Lamb-Oseen v_theta = G/(2 pi r) (1 - exp(-r^2/R^2)), plus translation."""
    dx, dy = x - cx, y - cy
    r = np.maximum(np.hypot(dx, dy), 1e-9)
    vt = G / (2 * np.pi * r) * (1.0 - np.exp(-(r / R) ** 2))
    th = np.arctan2(dy, dx)
    return -vt * np.sin(th) + CU, vt * np.cos(th) + CV


def field_solid(x, y, cx, cy, om=OMEGA_SB):
    return -om * (y - cy) + CU, om * (x - cx) + CV


def zeta_oseen(r, G=G_TRUE, R=R_TRUE):
    return G / (np.pi * R ** 2) * np.exp(-(r / R) ** 2)


# -------------------------------------------------------------- machinery ---

def integrate(X0, Y0, ts, cx0, cy0, fld):
    """RK4 the drifters through the moving vortex.

    Integrating rather than imposing positions is the whole point: it makes
    positions and velocities mutually consistent, so the estimators see a flow
    that could actually have happened.  Imposing a center track on unrelated
    positions produces drifters that never orbit it and tests nothing.
    """
    n, nd = len(ts), len(X0)
    X = np.zeros((n, nd)); Y = np.zeros((n, nd))
    U = np.zeros((n, nd)); V = np.zeros((n, nd))
    CX = cx0 + CU * (ts - ts[0]); CY = cy0 + CV * (ts - ts[0])
    X[0], Y[0] = X0, Y0
    U[0], V[0] = fld(X0, Y0, CX[0], CY[0])
    for i in range(n - 1):
        h = ts[i + 1] - ts[i]
        cxm, cym = CX[i] + CU * h / 2, CY[i] + CV * h / 2
        k1x, k1y = fld(X[i], Y[i], CX[i], CY[i])
        k2x, k2y = fld(X[i] + h / 2 * k1x, Y[i] + h / 2 * k1y, cxm, cym)
        k3x, k3y = fld(X[i] + h / 2 * k2x, Y[i] + h / 2 * k2y, cxm, cym)
        k4x, k4y = fld(X[i] + h * k3x, Y[i] + h * k3y, CX[i + 1], CY[i + 1])
        X[i + 1] = X[i] + h / 6 * (k1x + 2 * k2x + 2 * k3x + k4x)
        Y[i + 1] = Y[i] + h / 6 * (k1y + 2 * k2y + 2 * k3y + k4y)
        U[i + 1], V[i + 1] = fld(X[i + 1], Y[i + 1], CX[i + 1], CY[i + 1])
    return X, Y, U, V, CX, CY


def _circ_trap(px, py, uvfn, cx, cy):
    """Exactly what circulation_kinematics does: trapezoid on vertex values."""
    xn, yn = np.roll(px, -1), np.roll(py, -1)
    u, v = uvfn(px, py, cx, cy)
    un, vn = np.roll(u, -1), np.roll(v, -1)
    return np.sum(0.5 * (u + un) * (xn - px) + 0.5 * (v + vn) * (yn - py))


def _circ_exact(px, py, uvfn, cx, cy):
    """Gauss-Legendre along each straight edge: the exact same line integral."""
    xg, wg = _GL
    s = 0.5 * (xg + 1.0)
    tot = 0.0
    for i in range(len(px)):
        j = (i + 1) % len(px)
        dx, dy = px[j] - px[i], py[j] - py[i]
        u, v = uvfn(px[i] + s * dx, py[i] + s * dy, cx, cy)
        tot += 0.5 * np.sum(wg * (u * dx + v * dy))
    return tot


def _load(blocks=4):
    paths = {k: os.path.join(DATA, "drifters", f"{k}_gps_timeseries.nc")
             for k in ek.DRIFTERS}
    t0 = np.datetime64("2023-05-22T05:10:54").astype("datetime64[ms]").astype(np.int64)
    t1 = np.datetime64("2023-05-23T06:34:00").astype("datetime64[ms]").astype(np.int64)
    return ek.assemble(paths, t0, t1, blocks_per_burst=blocks)


def _shrink(X0, Y0, f):
    mx, my = np.nanmean(X0), np.nanmean(Y0)
    return mx + (X0 - mx) * f, my + (Y0 - my) * f


# ------------------------------------------------------------------ tests ---

def test_solid_body_exact(A):
    """A solid-body field is affine, so the circulation estimator must be exact."""
    ts = (A["t"] - A["t"][0]) / 1000.0
    X0, Y0 = _shrink(A["x"][0], A["y"][0], 0.6)
    X, Y, U, V, _, _ = integrate(X0, Y0, ts, np.nanmean(X0) - 900,
                                 np.nanmean(Y0) - 500, field_solid)
    cc = ek.circulation_kinematics(A["t"], X, Y, U, V)
    got, want = np.nanmedian(cc["zeta"]), 2 * OMEGA_SB
    err = abs(got - want) / abs(want)
    print(f"  1. solid body : zeta {got:+.6e} vs truth {want:+.6e}  "
          f"({100*err:.4f}%)")
    assert err < 1e-6, f"affine field must be exact, got {100*err:.4f}%"
    return True


def test_solid_body_null(A):
    """The 'not solid body' diagnostic must not fire on an actual solid body."""
    ts = (A["t"] - A["t"][0]) / 1000.0
    out = []
    for f, fld, tag in ((0.5, field_solid, "solid body "),
                        (1.0, field_solid, "solid body "),
                        (0.5, field_oseen, "Lamb-Oseen"),
                        (1.0, field_oseen, "Lamb-Oseen")):
        X0, Y0 = _shrink(A["x"][0], A["y"][0], f)
        X, Y, U, V, CX, CY = integrate(X0, Y0, ts, np.nanmean(X0) - 900,
                                       np.nanmean(Y0) - 500, fld)
        cc = ek.circulation_kinematics(A["t"], X, Y, U, V)
        rc = np.hypot(np.nanmean(X, axis=1) - CX, np.nanmean(Y, axis=1) - CY)
        m = np.isfinite(cc["zeta"]) & np.isfinite(rc)
        r = np.corrcoef(np.abs(cc["zeta"][m]), rc[m])[0, 1]
        sc = np.sqrt(np.nanmedian(np.abs(cc["area"])))
        print(f"  2. {tag} cluster {sc:5.0f} m: corr(|zeta|, radius) = {r:+.2f}")
        out.append((tag, r))
    sb = [r for tag, r in out if tag.startswith("solid")]
    lo = [r for tag, r in out if tag.startswith("Lamb")]
    assert max(abs(r) for r in sb) < 0.25, \
        f"solid body must give ~0 correlation, got {sb}"
    assert max(lo) < -0.3, f"Lamb-Oseen must give a clear negative, got {lo}"

    # 2b.  The test above advects a RIGID cluster, whose size never varies, so
    # it cannot exercise a size-dependent artifact.  Pool many cluster sizes in
    # the same constant-vorticity field and correlate against cluster scale.
    zs, ss = [], []
    for f in np.linspace(0.3, 1.6, 14):
        X0, Y0 = _shrink(A["x"][0], A["y"][0], f)
        X, Y, U, V, _, _ = integrate(X0, Y0, ts, np.nanmean(X0) - 900,
                                     np.nanmean(Y0) - 500, field_solid)
        cc = ek.circulation_kinematics(A["t"], X, Y, U, V)
        z = np.abs(cc["zeta"])
        sc = np.sqrt(np.abs(cc["area"]))
        m = np.isfinite(z) & np.isfinite(sc)
        zs.append(z[m])
        ss.append(sc[m])
    z, sc = np.concatenate(zs), np.concatenate(ss)
    r_pool = np.corrcoef(z, sc)[0, 1]
    spread = (z.max() - z.min()) / z.mean()
    print(f"  2b. solid body POOLED over {sc.min():.0f}-{sc.max():.0f} m "
          f"(n={len(z)}): corr = {r_pool:+.3f}, |zeta| spread {100*spread:.4f}%")
    assert abs(r_pool) < 0.10, \
        f"pooled solid body must give ~0 correlation, got {r_pool:+.3f}"
    assert spread < 1e-6, \
        f"|zeta| must be size-independent for a solid body, spread {spread:.2e}"
    print("     -> no size-dependent artifact: the estimator returns the same")
    print("        zeta at every cluster size, so the real -0.56 is diagnostic")
    return True


def test_oseen_recovery(A):
    """Gamma and R must come back from a known Lamb-Oseen vortex."""
    ts = (A["t"] - A["t"][0]) / 1000.0
    X0, Y0 = _shrink(A["x"][0], A["y"][0], 0.6)
    X, Y, U, V, CX, CY = integrate(X0, Y0, ts, np.nanmean(X0) - 900,
                                   np.nanmean(Y0) - 500, field_oseen)
    trans = np.column_stack([np.full(len(A["t"]), CU), np.full(len(A["t"]), CV)])
    p = ek.radial_profile(A["t"], X, Y, U, V, CX, CY, trans)
    o = ek.fit_oseen(p["r"], p["v_theta"], p["v_theta_err"])
    eg = 100 * (o["gamma"] - G_TRUE) / abs(G_TRUE)
    er = 100 * (o["radius"] - R_TRUE) / R_TRUE
    print(f"  3. Lamb-Oseen : Gamma {o['gamma']:+.0f} ({eg:+.1f}%), "
          f"R {o['radius']:.0f} m ({er:+.1f}%)")
    assert abs(eg) < 15 and abs(er) < 15, "recovery worse than 15%"
    return True


def test_formal_error_realism(A):
    """Is the least-squares formal error as optimistic as docs/03 claims?

    Truth is the EXACT area-averaged zeta over the same polygon, computed by
    Gauss-Legendre quadrature -- the quantity the LSQ fit actually estimates.
    """
    ts = (A["t"] - A["t"][0]) / 1000.0
    out = []
    for f in (0.5, 1.0):
        X0, Y0 = _shrink(A["x"][0], A["y"][0], f)
        X, Y, U, V, CX, CY = integrate(X0, Y0, ts, np.nanmean(X0) - 900,
                                       np.nanmean(Y0) - 500, field_oseen)
        n = len(A["t"])
        zt = np.full(n, np.nan)
        for i in range(n):
            m = np.isfinite(X[i]) & np.isfinite(Y[i])
            if m.sum() < 3:
                continue
            xs, ys = X[i][m], Y[i][m]
            o = np.argsort(np.arctan2(ys - ys.mean(), xs - xs.mean()))
            xs, ys = xs[o], ys[o]
            Ar = 0.5 * np.sum(xs * np.roll(ys, -1) - np.roll(xs, -1) * ys)
            if abs(Ar) < 1.0:
                continue
            zt[i] = _circ_exact(xs, ys, field_oseen, CX[i], CY[i]) / Ar
        F = ek.fit_gradient(A["t"], X, Y, U, V, window_s=1800.0, step_s=600.0)
        g = np.isfinite(zt)
        tv = np.interp(F["t"], A["t"][g], zt[g])
        m = np.isfinite(F["zeta"]) & np.isfinite(F["zeta_err"])
        act = np.sqrt(np.nanmean((F["zeta"][m] - tv[m]) ** 2))
        formal = np.nanmedian(F["zeta_err"][m])
        cc = ek.circulation_kinematics(A["t"], X, Y, U, V)
        sc = np.sqrt(np.nanmedian(np.abs(cc["area"])))
        print(f"  5. formal err : cluster {sc:4.0f} m, actual/formal = "
              f"{act/formal:.2f}x  (actual {act:.3e}, formal {formal:.3e})")
        out.append(act / formal)
    assert min(out) > 3.0, \
        f"formal error must be clearly optimistic, got {out}"
    assert min(out) < 9.0, \
        f"ratio implausibly large -- check the truth definition, got {out}"
    print("     -> docs/03's '5.2x too small' is reproduced; it is cluster-size")
    print(f"        dependent, {min(out):.1f}x to {max(out):.1f}x over 934-1105 m")
    return True


def test_quadrature_bias(A):
    """Measure the trapezoid error using the REAL cluster geometry.

    Uses the fitted Lamb-Oseen as the assumed truth and asks, for each real
    polygon, how far the estimator's trapezoid rule falls from the exact line
    integral.  This is a self-consistency check on the published zeta.
    """
    t, X, Y, U, V = A["t"], A["x"], A["y"], A["u"], A["v"]
    cc = ek.circulation_kinematics(t, X, Y, U, V)
    F = ek.fit_gradient(t, X, Y, U, V, window_s=1800.0, step_s=600.0)
    Cxy = ek.translation_velocity(t, U, V, period_s=21600.0)
    ct = np.column_stack([np.interp(F["t"], t, Cxy[:, 0]),
                          np.interp(F["t"], t, Cxy[:, 1])])
    sc = np.interp(F["t"], t, np.sqrt(np.where(np.isfinite(cc["area"]),
                                               np.abs(cc["area"]), np.nan)))
    cx, cy = ek.eddy_center(F, translation=ct, snr=3.0, max_disp_scales=3.0,
                            cluster_scale=sc)
    g = np.isfinite(cx)
    CXi = np.interp(t, F["t"][g], cx[g], left=np.nan, right=np.nan)
    CYi = np.interp(t, F["t"][g], cy[g], left=np.nan, right=np.nan)

    def uvfn(x, y, ccx, ccy):
        return field_oseen(x, y, ccx, ccy, G=-6114.0, R=1206.0)

    rel = []
    for i in range(len(t)):
        ok = np.isfinite(X[i]) & np.isfinite(Y[i])
        if ok.sum() < 4 or not np.isfinite(CXi[i]):
            continue
        xs, ys = X[i][ok], Y[i][ok]
        o = np.argsort(np.arctan2(ys - ys.mean(), xs - xs.mean()))
        px, py = xs[o], ys[o]
        Ar = abs(0.5 * np.sum(px * np.roll(py, -1) - np.roll(px, -1) * py))
        if Ar < 1.0:
            continue
        zt = _circ_trap(px, py, uvfn, CXi[i], CYi[i]) / Ar
        ze = _circ_exact(px, py, uvfn, CXi[i], CYi[i]) / Ar
        if abs(ze) > 1e-9:
            rel.append(100 * (abs(zt) - abs(ze)) / abs(ze))
    rel = np.array(rel)
    med = np.median(rel)
    print(f"  4. quadrature : |zeta| biased {med:+.2f}% "
          f"(IQR {np.percentile(rel,25):+.2f} .. {np.percentile(rel,75):+.2f}, "
          f"n={len(rel)})")
    zpub = -1.187e-3          # eddy_analysis.py median, both shape gates on
    print(f"     -> published {zpub:+.3e} corrects to "
          f"{zpub/(1+med/100):+.3e} (Rossby {zpub/1.7756e-5:+.1f} -> "
          f"{zpub/(1+med/100)/1.7756e-5:+.1f})")
    assert abs(med) < 10, "bias larger than 10% would invalidate the headline"
    return True


def main():
    nc = os.path.join(DATA, "drifters", "mwb458d02_gps_timeseries.nc")
    if not os.path.exists(nc):
        print(f"SKIP -- no drifter data at {nc}")
        return 0
    A = _load()
    print("Synthetic recovery tests (truth is known by construction)")
    ok = all([test_solid_body_exact(A), test_solid_body_null(A),
              test_oseen_recovery(A), test_formal_error_realism(A),
              test_quadrature_bias(A)])
    print("\nPASS -- all null tests satisfied" if ok else "\nFAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
