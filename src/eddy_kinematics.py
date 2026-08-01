#! /usr/bin/env python3

"""Velocity-gradient kinematics of the Peleliu tip vortex from four mwb drifters.

WHAT THIS DOES DIFFERENTLY FROM explore_drifter_paths.m
-------------------------------------------------------
explore_drifter_paths.m models each drifter's *position* as a shared eddy
centroid plus a per-drifter radial offset, taking the orbital direction from
the drifter's own smoothed velocity (so the center is assumed to lie +/-90 deg
off the velocity vector).  It fits 2 + nDrifters parameters per window with
fminsearch.  Vorticity is then only obtainable indirectly, via the fitted
radius and orbital speed, and the fit carries a sign degeneracy (the center may
sit on either side of the velocity vector, so radii can come out negative).

Here the target quantity -- vorticity -- is estimated directly from the
measured velocity field.  Over a cluster small compared with the eddy, the
velocity field is approximately affine:

    u(x,y,t) = U0 + Ut*(t-tc) + dudx*(x-xc) + dudy*(y-yc)
    v(x,y,t) = V0 + Vt*(t-tc) + dvdx*(x-xc) + dvdy*(y-yc)

Every drifter sample in a time window is one row of a linear system, so the
velocity-gradient tensor comes from ordinary least squares, in closed form,
with a full covariance matrix.  From the tensor:

    vorticity        zeta  = dvdx - dudy
    divergence       delta = dudx + dvdy
    normal strain    sig_n = dudx - dvdy
    shear strain     sig_s = dvdx + dudy
    Okubo-Weiss      OW    = sig_n^2 + sig_s^2 - zeta^2   (OW < 0 => rotation)

Advantages over the geometric fit: vorticity is measured rather than inferred;
divergence and strain come free; the system is linear so there is no starting
guess, no fminsearch tolerance, and no sign degeneracy; and formal uncertainty
follows from the covariance.  The eddy center is recovered from the same fit as
the elliptic critical point (see eddy_center), rather than being the primary
fitted parameter.

INDEPENDENT CROSS-CHECK
-----------------------
Vorticity and divergence are also computed geometrically, by Stokes' and the
divergence theorem around the drifter quadrilateral:

    zeta_circ  = (1/A) * closed_integral(u . dl)
    delta_circ = (1/A) * closed_integral(u . n dl)

For a strictly affine field these agree with the least-squares values exactly,
so any systematic difference measures how badly the affine assumption fails
across the cluster.  Divergence gets a third, fully independent estimate from
the cluster area itself, delta_area = d(ln A)/dt, which uses only positions and
no velocities at all.  Agreement among three estimators built on different
information is the main evidence that the numbers mean something.

PREPROCESSING
-------------
Raw 2 Hz velocity is dominated by surface wave orbital motion (spectral peak
near 4.4 s, essentially all wave energy below 30 s period).  Each 2048-sample
burst is therefore block-averaged: bursts are internally gap-free, so averaging
within a burst never straddles the 35-202 s inter-burst gaps.  With the default
4 blocks per burst the output cadence is ~256 s, far faster than the eddy
orbital period (hours) and far slower than the waves.

USAGE
-----
    python3 eddy_kinematics.py --input . --output ./eddy_out \\
        --start 2023-05-22T05:10:54 --end 2023-05-23T06:34:00
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2 - WGS84_F)
OMEGA = 7.292115e-5           # Earth rotation rate, rad/s

DRIFTERS = ("mwb458d02", "mwb788d01", "mwb790d01", "mwb793d02")


# ------------------------------------------------------------------ frame ---

def local_frame(lat, lon, lat0, lon0):
    """Geodetic -> local east/north tangent plane (meters) about (lat0, lon0).

    Uses the WGS-84 meridional and prime-vertical radii of curvature at lat0.
    Exact to better than a millimeter over the few-km extent of this cluster.
    """
    p = np.radians(lat0)
    s = np.sin(p)
    rm = WGS84_A * (1 - WGS84_E2) / (1 - WGS84_E2 * s * s) ** 1.5   # north
    rn = WGS84_A / np.sqrt(1 - WGS84_E2 * s * s)                    # east
    x = np.radians(lon - lon0) * rn * np.cos(p)
    y = np.radians(lat - lat0) * rm
    return x, y


def coriolis(lat0):
    return 2.0 * OMEGA * np.sin(np.radians(lat0))


# ------------------------------------------------------------------- data ---

def load_blocks(path, t0_ms=None, t1_ms=None, blocks_per_burst=4):
    """Read one drifter's netCDF and block-average each burst.

    Averaging happens strictly inside a burst, which is gap-free, so no average
    ever spans an inter-burst gap.  Returns a dict of 1-D arrays ordered in
    time: t (ms), lat, lon, u, v, plus n (samples per block).
    """
    import netCDF4
    ds = netCDF4.Dataset(path)
    for v in ds.variables.values():
        v.set_auto_mask(False)
    T = ds.variables["time"][:].astype(np.int64)
    nb, ns = T.shape
    if ns % blocks_per_burst:
        raise ValueError(f"{ns} samples per burst is not divisible by "
                         f"{blocks_per_burst}")
    keep = np.ones(nb, dtype=bool)
    if t0_ms is not None:
        keep &= T[:, -1] >= t0_ms
    if t1_ms is not None:
        keep &= T[:, 0] <= t1_ms
    if not keep.any():
        raise ValueError(f"{path}: no bursts inside the requested window")

    out = {}
    m = ns // blocks_per_burst
    for name, dtype in (("time", np.float64), ("lat", np.float64),
                        ("lon", np.float64), ("u", np.float64),
                        ("v", np.float64)):
        a = np.asarray(ds.variables[name][:], dtype)[keep]
        out[name] = a.reshape(-1, blocks_per_burst, m).mean(axis=2).ravel()
    ds.close()
    out["t"] = out.pop("time")
    out["n"] = m
    # keep only blocks whose center falls in the window
    sel = np.ones(len(out["t"]), dtype=bool)
    if t0_ms is not None:
        sel &= out["t"] >= t0_ms
    if t1_ms is not None:
        sel &= out["t"] <= t1_ms
    for k in ("t", "lat", "lon", "u", "v"):
        out[k] = out[k][sel]
    return out


def assemble(paths, t0_ms=None, t1_ms=None, blocks_per_burst=4, max_gap_s=1800.0):
    """Load all drifters, build a common local frame and a common time grid.

    Each drifter's block series is linearly interpolated onto the common grid.
    Grid points further than ``max_gap_s`` from a real observation of a given
    drifter are marked invalid for that drifter rather than extrapolated.
    """
    d = {k: load_blocks(p, t0_ms, t1_ms, blocks_per_burst) for k, p in paths.items()}
    names = list(d)
    lat0 = float(np.mean([np.mean(v["lat"]) for v in d.values()]))
    lon0 = float(np.mean([np.mean(v["lon"]) for v in d.values()]))
    for v in d.values():
        v["x"], v["y"] = local_frame(v["lat"], v["lon"], lat0, lon0)

    lo = max(v["t"].min() for v in d.values())
    hi = min(v["t"].max() for v in d.values())
    step = np.median(np.concatenate([np.diff(v["t"]) for v in d.values()]))
    grid = np.arange(lo, hi + 1, step)

    nT, nD = len(grid), len(names)
    X = np.full((nT, nD), np.nan)
    Y, U, V = X.copy(), X.copy(), X.copy()
    for j, k in enumerate(names):
        v = d[k]
        for dst, src in ((X, "x"), (Y, "y"), (U, "u"), (V, "v")):
            dst[:, j] = np.interp(grid, v["t"], v[src])
        near = np.abs(grid[:, None] - v["t"][None, :]).min(axis=1) / 1000.0
        bad = near > max_gap_s
        for dst in (X, Y, U, V):
            dst[bad, j] = np.nan
    return dict(names=names, t=grid, x=X, y=Y, u=U, v=V,
                lat0=lat0, lon0=lon0, f=coriolis(lat0), raw=d, step_ms=step)


# ------------------------------------------------- least-squares gradient ---

def fit_gradient(t, x, y, u, v, window_s=3600.0, step_s=None, time_term=True,
                 min_drifters=3):
    """Sliding-window least-squares fit of the affine velocity field.

    The u and v equations share a design matrix but no parameters, so each is
    an independent linear system; they are solved together for convenience.

    Returns a dict of time series: zeta, delta, sig_n, sig_s, OW, their formal
    1-sigma uncertainties, the window-mean flow, condition number and residual.
    """
    step_s = step_s or window_s / 8.0
    tc = np.arange(t.min(), t.max() + 1, step_s * 1000.0)
    half = window_s * 1000.0 / 2.0
    ncol = 4 if time_term else 3

    keys = ("zeta", "delta", "sig_n", "sig_s", "OW", "U0", "V0",
            "zeta_err", "delta_err", "sig_n_err", "sig_s_err",
            "cond", "rms", "ndrift", "nobs", "xc", "yc")
    out = {k: np.full(len(tc), np.nan) for k in keys}
    grad = np.full((len(tc), 2, 2), np.nan)
    gcov = np.full((len(tc), 6, 6), np.nan)

    for i, t0 in enumerate(tc):
        m = np.abs(t - t0) <= half
        if not m.any():
            continue
        tt = (t[m] - t0) / 1000.0
        xx, yy, uu, vv = x[m], y[m], u[m], v[m]
        ok = np.isfinite(xx) & np.isfinite(yy) & np.isfinite(uu) & np.isfinite(vv)
        ndrift = int((ok.any(axis=0)).sum())
        if ndrift < min_drifters or ok.sum() < ncol + 2:
            continue
        T = np.repeat(tt[:, None], xx.shape[1], axis=1)[ok]
        xs, ys = xx[ok], yy[ok]
        xbar, ybar = xs.mean(), ys.mean()
        dx, dy = xs - xbar, ys - ybar
        cols = [np.ones_like(dx), dx, dy] + ([T] if time_term else [])
        M = np.column_stack(cols)
        # Columns carry different units (1, m, m, s), so an unscaled design
        # matrix reports a condition number dominated by that unit mismatch
        # rather than by cluster geometry.  Scale to unit norm, solve, then
        # undo the scaling; this leaves the fit unchanged but makes cond(M) a
        # meaningful statement about how well the drifters span the plane.
        sc = np.linalg.norm(M, axis=0)
        sc[sc == 0] = 1.0
        Ms = M / sc

        # solve both components against the same design matrix
        B = np.column_stack([uu[ok], vv[ok]])
        coefs, *_ = np.linalg.lstsq(Ms, B, rcond=None)
        coef = coefs / sc[:, None]
        res = B - M @ coef
        dof = max(M.shape[0] - ncol, 1)
        s2 = (res ** 2).sum(axis=0) / dof                # per component
        MtM_inv = np.linalg.pinv(Ms.T @ Ms) / np.outer(sc, sc)

        U0, dudx, dudy = coef[0, 0], coef[1, 0], coef[2, 0]
        V0, dvdx, dvdy = coef[0, 1], coef[1, 1], coef[2, 1]
        vdudx, vdudy = s2[0] * MtM_inv[1, 1], s2[0] * MtM_inv[2, 2]
        vdvdx, vdvdy = s2[1] * MtM_inv[1, 1], s2[1] * MtM_inv[2, 2]

        zeta, delta = dvdx - dudy, dudx + dvdy
        sn, ss = dudx - dvdy, dvdx + dudy
        out["zeta"][i], out["delta"][i] = zeta, delta
        out["sig_n"][i], out["sig_s"][i] = sn, ss
        out["OW"][i] = sn ** 2 + ss ** 2 - zeta ** 2
        # u and v systems are independent, so their variances simply add
        out["zeta_err"][i] = np.sqrt(vdvdx + vdudy)
        out["delta_err"][i] = np.sqrt(vdudx + vdvdy)
        out["sig_n_err"][i] = np.sqrt(vdudx + vdvdy)
        out["sig_s_err"][i] = np.sqrt(vdvdx + vdudy)
        out["U0"][i], out["V0"][i] = U0, V0
        out["cond"][i] = np.linalg.cond(Ms)
        out["rms"][i] = np.sqrt((res ** 2).mean())
        out["ndrift"][i], out["nobs"][i] = ndrift, ok.sum()
        out["xc"][i], out["yc"][i] = xbar, ybar
        grad[i] = [[dudx, dudy], [dvdx, dvdy]]

    out["t"] = tc
    out["grad"] = grad
    return out


def eddy_center(fit, translation=None, snr=3.0, max_disp_scales=3.0,
                cluster_scale=None):
    """Locate the elliptic critical point of each windowed affine fit.

    In a frame moving with the eddy at velocity c, the center is where the
    velocity vanishes:  A (x - xbar) = c - U0, so x = xbar + A^-1 (c - U0).
    The point is a genuine vortex center (rather than a saddle) only when A has
    complex eigenvalues, equivalently when Okubo-Weiss is negative; elsewhere
    the result is returned as NaN.
    """
    n = len(fit["t"])
    cx = np.full(n, np.nan)
    cy = np.full(n, np.nan)
    if translation is None:
        translation = np.zeros((n, 2))
    for i in range(n):
        A = fit["grad"][i]
        if not np.all(np.isfinite(A)) or fit["OW"][i] >= 0:
            continue
        # A^-1 scales as 1/|zeta|, so where the rotation is weak or poorly
        # determined the inversion amplifies a small velocity mismatch into a
        # huge displacement.  Require the vorticity to be resolved before
        # inverting at all.
        if not np.isfinite(fit["zeta_err"][i]) or \
                abs(fit["zeta"][i]) < snr * fit["zeta_err"][i]:
            continue
        rhs = np.array([translation[i, 0] - fit["U0"][i],
                        translation[i, 1] - fit["V0"][i]])
        try:
            if abs(np.linalg.det(A)) < 1e-14:
                continue
            d = np.linalg.solve(A, rhs)
        except np.linalg.LinAlgError:
            continue
        if cluster_scale is not None and np.isfinite(cluster_scale[i]):
            # a center many cluster-widths away is not constrained by drifters
            # that all sit on one side of it
            if np.hypot(*d) > max_disp_scales * cluster_scale[i]:
                continue
        cx[i], cy[i] = fit["xc"][i] + d[0], fit["yc"][i] + d[1]
    return cx, cy


# ------------------------------------------------- geometric cross-checks ---

def _order_ccw(x, y):
    a = np.arctan2(y - y.mean(), x - x.mean())
    return np.argsort(a)


def polygon_quality(x, y):
    """Isoperimetric quotient 4*pi*A / P**2 of a polygon: 1 for a circle,
    ~0.785 for a square, -> 0 as the vertices approach collinearity.

    Vorticity from circulation is Gamma/A, so a sliver-shaped cluster divides a
    small, noisy circulation by a vanishing area and the estimate diverges.
    This is the natural quantity to gate on: it is scale-free, so it rejects
    degenerate *shapes* without rejecting a merely small cluster.
    """
    n = len(x)
    A = abs(0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))
    P = np.sum(np.hypot(np.roll(x, -1) - x, np.roll(y, -1) - y))
    return 4 * np.pi * A / P ** 2 if P > 0 else 0.0


def circulation_kinematics(t, x, y, u, v, min_quality=0.10):
    """Vorticity and divergence from contour integrals around the drifters.

    Stokes' theorem gives the area-averaged vorticity as the circulation
    divided by the enclosed area; the divergence theorem gives the
    area-averaged divergence as the outward flux divided by the area.  Both use
    the drifter polygon directly and share no machinery with the least-squares
    fit, so they are an independent check on it.  Requires all four drifters.
    """
    n = len(t)
    zeta = np.full(n, np.nan)
    delta = np.full(n, np.nan)
    area = np.full(n, np.nan)
    aspect = np.full(n, np.nan)
    qual = np.full(n, np.nan)
    for i in range(n):
        ok = (np.isfinite(x[i]) & np.isfinite(y[i])
              & np.isfinite(u[i]) & np.isfinite(v[i]))
        if ok.sum() < 3:
            continue
        xs, ys, us, vs = x[i][ok], y[i][ok], u[i][ok], v[i][ok]
        o = _order_ccw(xs, ys)
        xs, ys, us, vs = xs[o], ys[o], us[o], vs[o]
        xn, yn = np.roll(xs, -1), np.roll(ys, -1)
        un, vn = np.roll(us, -1), np.roll(vs, -1)
        A = 0.5 * np.sum(xs * yn - xn * ys)            # shoelace, signed
        qual[i] = polygon_quality(xs, ys)
        if abs(A) < 1.0 or qual[i] < min_quality:
            continue                                   # degenerate cluster
        dx, dy = xn - xs, yn - ys
        circ = np.sum(0.5 * (us + un) * dx + 0.5 * (vs + vn) * dy)
        flux = np.sum(0.5 * (us + un) * dy - 0.5 * (vs + vn) * dx)
        zeta[i], delta[i], area[i] = circ / A, flux / A, A
        # cluster shape: ratio of principal axes of the position covariance
        c = np.cov(np.vstack([xs, ys]))
        w = np.linalg.eigvalsh(c)
        aspect[i] = np.sqrt(max(w, default=np.nan) / max(w.min(), 1e-12)) \
            if np.ndim(w) == 0 else np.sqrt(w.max() / max(w.min(), 1e-12))
    return dict(zeta=zeta, delta=delta, area=area, aspect=aspect, quality=qual)


def circulation_jackknife(t, x, y, u, v):
    """Leave-one-drifter-out spread of the circulation estimates.

    Four drifters give four distinct triangles.  For a strictly affine field
    every triangle returns the same vorticity, so the spread across them
    measures the combined effect of measurement noise and of curvature in the
    velocity field across the cluster.  This is an empirical uncertainty that
    assumes neither linearity nor a noise model, unlike the least-squares
    formal error.
    """
    n, nd = x.shape
    zj = np.full((n, nd), np.nan)
    dj = np.full((n, nd), np.nan)
    for k in range(nd):
        keep = [j for j in range(nd) if j != k]
        c = circulation_kinematics(t, x[:, keep], y[:, keep],
                                   u[:, keep], v[:, keep])
        zj[:, k], dj[:, k] = c["zeta"], c["delta"]
    def rsd(a):
        med = np.nanmedian(a, axis=1, keepdims=True)
        return 1.4826 * np.nanmedian(np.abs(a - med), axis=1)

    # median and MAD, not mean and std: with only four triangles a single
    # near-degenerate one would otherwise dominate both.
    return dict(zeta=zj, delta=dj,
                zeta_spread=rsd(zj), delta_spread=rsd(dj),
                zeta_mean=np.nanmedian(zj, axis=1),
                n_ok=np.isfinite(zj).sum(axis=1))


def rotation_rate(t, x, y):
    """Cumulative rotation of the drifter constellation, positions only.

    Tracks the mean angle of the drifters about their own centroid, unwrapped,
    so the total turning over the record can be compared with the time integral
    of vorticity: for solid-body rotation d(theta)/dt = zeta/2.  Uses no
    velocity data at all, so it is independent of every other estimator here.
    """
    n, nd = x.shape
    xc, yc = np.nanmean(x, axis=1), np.nanmean(y, axis=1)
    ang = np.arctan2(y - yc[:, None], x - xc[:, None])
    th = np.full(n, np.nan)
    for j in range(nd):
        a = np.unwrap(ang[:, j])
        th = a if j == 0 else th + a
    th /= nd
    ts = (t - t[0]) / 1000.0
    ok = np.isfinite(th)
    om = np.full(n, np.nan)
    om[ok] = np.gradient(th[ok], ts[ok])
    return dict(theta=th, omega=om, turns=(th[ok][-1] - th[ok][0]) / (2 * np.pi))


def divergence_from_area(t, area, smooth_s=3600.0):
    """delta = d(ln A)/dt -- a positions-only divergence estimate.

    Shares no velocity information with either of the other two estimators.
    """
    ok = np.isfinite(area) & (area > 0)
    out = np.full(len(t), np.nan)
    if ok.sum() < 5:
        return out
    ts = t[ok] / 1000.0
    la = np.log(area[ok])
    w = max(int(smooth_s / np.median(np.diff(ts))) | 1, 3)
    k = np.ones(w) / w
    lasm = np.convolve(np.pad(la, w // 2, mode="edge"), k, "valid")[:len(la)]
    out[ok] = np.gradient(lasm, ts)
    return out


def translation_velocity(t, u, v, period_s=21600.0):
    """Eddy translation, as the drifter-mean velocity low-passed over ``period_s``.

    At any single instant, a uniform background flow and a displacement of the
    vortex center are exactly degenerate: omega*zhat x (x - c) contains the
    constant -omega*zhat x c, which is indistinguishable from translation.  So
    the center cannot be recovered from an instantaneous velocity field alone.
    Averaging the drifter velocities over several orbital periods cancels the
    rotational part and leaves the translation, which breaks the degeneracy.
    Every center estimate below is conditional on this assumption.
    """
    m = np.nanmean(u, axis=1), np.nanmean(v, axis=1)
    ts = (t - t[0]) / 1000.0
    w = max(int(period_s / np.median(np.diff(ts))) | 1, 3)
    k = np.ones(w) / w
    out = []
    for a in m:
        b = np.where(np.isfinite(a), a, np.nanmean(a))
        out.append(np.convolve(np.pad(b, w // 2, mode="edge"), k, "valid")[:len(b)])
    return np.column_stack(out)


# --------------------------------------------------------- radial profile ---

def radial_profile(t, x, y, u, v, cx, cy, trans, nbin=14, rmax=None):
    """Azimuthal velocity and local vorticity as functions of radius.

    Each circulation estimate above is a mean vorticity over whatever area the
    drifter polygon happened to span, so it mixes radii.  Referred to the
    fitted center and the eddy translation, every drifter sample instead
    becomes one (r, v_theta) pair, and pooling them over the record resolves
    the structure the area-average was hiding:

        Gamma(r) = 2 pi r v_theta(r)          circulation within radius r
        zeta(r)  = (1/r) d(r v_theta)/dr      local vorticity

    Radial velocity v_r is returned too; for a coherent, non-dispersing vortex
    it should be small compared with v_theta, which is a check rather than an
    assumption.
    """
    dx = x - cx[:, None]
    dy = y - cy[:, None]
    r = np.hypot(dx, dy)
    th = np.arctan2(dy, dx)
    du = u - trans[:, 0][:, None]
    dv = v - trans[:, 1][:, None]
    vt = -du * np.sin(th) + dv * np.cos(th)
    vr = du * np.cos(th) + dv * np.sin(th)

    ok = np.isfinite(r) & np.isfinite(vt)
    rr, vv, ww = r[ok], vt[ok], vr[ok]
    rmax = rmax or np.nanpercentile(rr, 97)
    edges = np.linspace(0, rmax, nbin + 1)
    mid = 0.5 * (edges[1:] + edges[:-1])
    vth = np.full(nbin, np.nan)
    vthe = np.full(nbin, np.nan)
    vrad = np.full(nbin, np.nan)
    cnt = np.zeros(nbin, int)
    for i in range(nbin):
        m = (rr >= edges[i]) & (rr < edges[i + 1])
        cnt[i] = m.sum()
        if cnt[i] < 20:
            continue
        vth[i] = np.median(vv[m])
        vthe[i] = 1.4826 * np.median(np.abs(vv[m] - vth[i])) / np.sqrt(cnt[i])
        vrad[i] = np.median(ww[m])
    good = np.isfinite(vth)
    zeta = np.full(nbin, np.nan)
    if good.sum() > 2:
        rv = mid * vth
        zeta[good] = np.gradient(rv[good], mid[good]) / mid[good]
    return dict(r=mid, v_theta=vth, v_theta_err=vthe, v_r=vrad,
                zeta=zeta, n=cnt, r_all=rr, vt_all=vv)


def fit_oseen(r, vth, w=None):
    """Fit a Lamb-Oseen vortex  v_theta = Gamma/(2 pi r) * (1 - exp(-r^2/R^2)).

    Returns circulation Gamma (m^2/s), core radius R (m) and the peak core
    vorticity Gamma/(pi R^2).  Brute-force over R with Gamma solved linearly at
    each R, so there is no starting guess and no local-minimum risk.
    """
    ok = np.isfinite(r) & np.isfinite(vth) & (r > 0)
    if ok.sum() < 3:
        return None
    rr, vv = r[ok], vth[ok]
    ww = np.ones_like(rr) if w is None else 1.0 / np.maximum(w[ok], 1e-9) ** 2
    best = None
    for R in np.linspace(0.05 * rr.max(), 2.0 * rr.max(), 400):
        basis = (1.0 - np.exp(-(rr / R) ** 2)) / (2 * np.pi * rr)
        g = np.sum(ww * basis * vv) / max(np.sum(ww * basis ** 2), 1e-30)
        resid = np.sum(ww * (vv - g * basis) ** 2)
        if best is None or resid < best[0]:
            best = (resid, g, R)
    _, G, R = best
    return dict(gamma=G, radius=R, zeta_core=G / (np.pi * R ** 2),
                vmax=G / (2 * np.pi * R) * (1 - np.exp(-1.0)) / 1.0)
