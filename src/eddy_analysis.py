#! /usr/bin/env python3

"""Run the drifter-cluster kinematics and write figures and a results file.

    python3 eddy_analysis.py --input . --output eddy_out \\
        --start 2023-05-22T05:10:54 --end 2023-05-23T06:34:00

Primary product is the vorticity time series of the Peleliu tip vortex, with
three estimators built on different information so they can be cross-checked,
and a leave-one-drifter-out uncertainty that assumes neither linearity of the
velocity field nor a noise model.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

import eddy_kinematics as ek
from paths import DATA


def _dt64(a):
    return np.asarray(a, dtype=np.int64).astype("datetime64[ms]")


def run(indir, outdir, t0, t1, blocks=4, window_s=1800.0, step_s=600.0,
        trans_period_s=21600.0):
    os.makedirs(outdir, exist_ok=True)
    paths = {k: os.path.join(indir, f"{k}_gps_timeseries.nc") for k in ek.DRIFTERS}
    D = ek.assemble(paths, t0, t1, blocks_per_burst=blocks)
    t, X, Y, U, V = D["t"], D["x"], D["y"], D["u"], D["v"]

    cc = ek.circulation_kinematics(t, X, Y, U, V)
    jk = ek.circulation_jackknife(t, X, Y, U, V)
    rot = ek.rotation_rate(t, X, Y)
    dA = ek.divergence_from_area(t, cc["area"], smooth_s=3600.0)
    F = ek.fit_gradient(t, X, Y, U, V, window_s=window_s, step_s=step_s)
    Cxy = ek.translation_velocity(t, U, V, period_s=trans_period_s)
    ct = np.column_stack([np.interp(F["t"], t, Cxy[:, 0]),
                          np.interp(F["t"], t, Cxy[:, 1])])
    scale = np.interp(F["t"], t, np.sqrt(np.where(np.isfinite(cc["area"]),
                                                  np.abs(cc["area"]), np.nan)))
    cx, cy = ek.eddy_center(F, translation=ct, snr=3.0, max_disp_scales=3.0,
                            cluster_scale=scale)
    print(f"eddy center: {np.isfinite(cx).sum()}/{len(cx)} windows pass the "
          f"rotation-resolved and displacement gates "
          f"({100*np.isfinite(cx).mean():.0f}%)")

    cxi = np.interp(t, F["t"][np.isfinite(cx)], cx[np.isfinite(cx)],
                    left=np.nan, right=np.nan)
    cyi = np.interp(t, F["t"][np.isfinite(cy)], cy[np.isfinite(cy)],
                    left=np.nan, right=np.nan)
    prof = ek.radial_profile(t, X, Y, U, V, cxi, cyi, Cxy)
    oseen = ek.fit_oseen(prof["r"], prof["v_theta"], prof["v_theta_err"])

    res = dict(t=t, x=X, y=Y, u=U, v=V, names=D["names"], f=D["f"],
               prof_r=prof["r"], prof_vt=prof["v_theta"],
               prof_vte=prof["v_theta_err"], prof_vr=prof["v_r"],
               prof_zeta=prof["zeta"], prof_n=prof["n"],
               oseen_gamma=(oseen or {}).get("gamma", np.nan),
               oseen_radius=(oseen or {}).get("radius", np.nan),
               oseen_zeta_core=(oseen or {}).get("zeta_core", np.nan),
               lat0=D["lat0"], lon0=D["lon0"],
               zeta_circ=cc["zeta"], delta_circ=cc["delta"],
               area=cc["area"], aspect=cc["aspect"], quality=cc["quality"],
               zeta_jk=jk["zeta_mean"], zeta_jk_spread=jk["zeta_spread"],
               delta_jk_spread=jk["delta_spread"],
               omega_pos=rot["omega"], theta_pos=rot["theta"], turns=rot["turns"],
               delta_area=dA, trans=Cxy,
               fit_t=F["t"], zeta_lsq=F["zeta"], zeta_lsq_err=F["zeta_err"],
               delta_lsq=F["delta"], sig_n=F["sig_n"], sig_s=F["sig_s"],
               OW=F["OW"], cond=F["cond"], rms=F["rms"], nobs=F["nobs"],
               center_x=cx, center_y=cy, window_s=window_s)
    np.savez_compressed(os.path.join(outdir, "eddy_kinematics.npz"),
                        **{k: v for k, v in res.items() if k != "names"})
    return res


# ------------------------------------------------------------------ plots ---

def figures(r, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"figure.dpi": 130, "font.size": 8,
                         "axes.grid": True, "grid.alpha": 0.25})
    T = _dt64(r["t"])
    FT = _dt64(r["fit_t"])
    f = r["f"]
    C = ["#3b6ea5", "#c1553b", "#4e8c58", "#8a6bab"]

    # 1 -- vorticity, the primary product
    fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True,
                           gridspec_kw={"height_ratios": [3, 1]})
    z, s = r["zeta_jk"], r["zeta_jk_spread"]
    ax[0].fill_between(T, z - s, z + s, color=C[0], alpha=.22,
                       label="leave-one-out spread (1$\\sigma$)")
    ax[0].plot(T, r["zeta_circ"], color=C[0], lw=1.1,
               label="circulation / Stokes (primary)")
    ax[0].plot(FT, r["zeta_lsq"], color=C[1], lw=1.1, alpha=.9,
               label=f"least-squares affine fit ({r['window_s']/3600:.1f} h window)")
    ax[0].plot(T, 2 * r["omega_pos"], color=C[2], lw=.9, alpha=.8,
               label="2 x constellation rotation (positions only)")
    ax[0].axhline(0, color="k", lw=.6)
    ax[0].axhline(f, color="0.4", ls=":", lw=.8)
    ax[0].text(T[3], f, "  +f", va="bottom", fontsize=7, color="0.4")
    ax[0].set_ylabel("vorticity $\\zeta$  (s$^{-1}$)")
    ax[0].legend(loc="upper right", fontsize=7, framealpha=.9)
    ax[0].set_title("Peleliu tip vortex: relative vorticity from four drifters")
    sec = ax[0].secondary_yaxis("right", functions=(lambda a: a / f, lambda a: a * f))
    sec.set_ylabel("Rossby number  $\\zeta/f$")
    ax[1].plot(T, r["area"] / 1e6, color="0.35", lw=1)
    ax[1].set_ylabel("cluster area\n(km$^2$)")
    ax[1].set_xlabel("2023 UTC")
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "vorticity.png"))
    plt.close(fig)

    # 2 -- full kinematic decomposition
    fig, ax = plt.subplots(3, 1, figsize=(9, 7), sharex=True)
    ax[0].plot(T, r["delta_circ"], color=C[0], lw=1, label="contour flux")
    ax[0].plot(T, r["delta_area"], color=C[2], lw=1, label="d(ln A)/dt")
    ax[0].plot(FT, r["delta_lsq"], color=C[1], lw=1, label="least squares")
    ax[0].axhline(0, color="k", lw=.6)
    ax[0].set_ylabel("divergence $\\delta$ (s$^{-1}$)")
    ax[0].legend(fontsize=7, ncol=3)
    ax[0].set_title("Divergence, strain and Okubo-Weiss")
    ax[1].plot(FT, r["sig_n"], color=C[0], lw=1, label="normal $\\sigma_n$")
    ax[1].plot(FT, r["sig_s"], color=C[1], lw=1, label="shear $\\sigma_s$")
    ax[1].plot(FT, np.hypot(r["sig_n"], r["sig_s"]), color="k", lw=1,
               label="total strain")
    ax[1].plot(FT, np.abs(r["zeta_lsq"]), color=C[2], lw=1, ls="--",
               label="|$\\zeta$|")
    ax[1].set_ylabel("strain (s$^{-1}$)")
    ax[1].legend(fontsize=7, ncol=4)
    ax[2].plot(FT, r["OW"], color="k", lw=1)
    ax[2].axhline(0, color="r", lw=.8)
    ax[2].fill_between(FT, r["OW"], 0, where=r["OW"] < 0, color=C[0], alpha=.25)
    ax[2].set_ylabel("Okubo-Weiss\n(s$^{-2}$)")
    ax[2].set_xlabel("2023 UTC")
    ax[2].text(.01, .08, "OW < 0: rotation dominates strain", transform=ax[2].transAxes,
               fontsize=7, color=C[0])
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "kinematics.png"))
    plt.close(fig)

    # 3 -- tracks and the inferred center
    fig, ax = plt.subplots(1, 2, figsize=(10, 5))
    for j, n in enumerate(r["names"]):
        ax[0].plot(r["x"][:, j] / 1e3, r["y"][:, j] / 1e3, lw=.8, color=C[j], label=n)
    ok = np.isfinite(r["center_x"])
    ax[0].plot(r["center_x"][ok] / 1e3, r["center_y"][ok] / 1e3, "k.-", ms=2, lw=.8,
               label="eddy center")
    ax[0].set_aspect("equal")
    ax[0].set_xlabel("east (km)"); ax[0].set_ylabel("north (km)")
    ax[0].legend(fontsize=7); ax[0].set_title("Earth frame")
    cxi = np.interp(r["t"], r["fit_t"][ok], r["center_x"][ok])
    cyi = np.interp(r["t"], r["fit_t"][ok], r["center_y"][ok])
    for j, n in enumerate(r["names"]):
        ax[1].plot((r["x"][:, j] - cxi), (r["y"][:, j] - cyi), lw=.7, color=C[j])
    ax[1].plot(0, 0, "k+", ms=12)
    ax[1].set_aspect("equal")
    ax[1].set_xlabel("east of center (m)"); ax[1].set_ylabel("north of center (m)")
    ax[1].set_title("Eddy frame")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "tracks_center.png"))
    plt.close(fig)

    # 4 -- estimator agreement and fit quality
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.4))
    zc = np.interp(r["fit_t"], r["t"], r["zeta_circ"])
    m = np.isfinite(zc) & np.isfinite(r["zeta_lsq"])
    ax[0].scatter(zc[m], r["zeta_lsq"][m], s=6, color=C[0], alpha=.6)
    lim = [min(zc[m].min(), r["zeta_lsq"][m].min()),
           max(zc[m].max(), r["zeta_lsq"][m].max())]
    ax[0].plot(lim, lim, "k--", lw=.8)
    ax[0].set_xlabel("$\\zeta$ circulation"); ax[0].set_ylabel("$\\zeta$ least squares")
    ax[0].set_title(f"r = {np.corrcoef(zc[m], r['zeta_lsq'][m])[0,1]:.2f}")
    ax[1].plot(T, r["quality"], color="0.3", lw=.9)
    ax[1].set_ylabel("polygon quality  4$\\pi$A/P$^2$")
    ax[1].axhline(0.10, color="r", ls=":", lw=.8)
    ax[1].set_title("cluster shape (<0.10 rejected)")
    rel = np.abs(r["zeta_jk_spread"] / r["zeta_jk"])
    ax[2].plot(T, 100 * rel, color=C[1], lw=.9)
    ax[2].set_ylabel("jackknife spread (% of $|\\zeta|$)")
    ax[2].set_ylim(0, 200)
    ax[2].set_title("uncertainty in $\\zeta$")
    for a in (ax[1], ax[2]):
        a.set_xlabel("2023 UTC")
        for lb in a.get_xticklabels():
            lb.set_rotation(30); lb.set_ha("right")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "diagnostics.png"))
    plt.close(fig)

    # 5 -- radial structure of the vortex
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))
    rr, vt, ve = r["prof_r"], r["prof_vt"], r["prof_vte"]
    ax[0].errorbar(rr, vt, yerr=ve, fmt="o-", ms=3, lw=1, color=C[0],
                   label="$v_\\theta$ observed")
    ax[0].plot(rr, r["prof_vr"], "s--", ms=3, lw=.8, color="0.55",
               label="$v_r$ (should be ~0)")
    if np.isfinite(r["oseen_radius"]):
        G, R = r["oseen_gamma"], r["oseen_radius"]
        rf = np.linspace(rr[np.isfinite(vt)].min(), rr.max(), 200)
        ax[0].plot(rf, G / (2 * np.pi * rf) * (1 - np.exp(-(rf / R) ** 2)),
                   "-", lw=1.2, color=C[1],
                   label=f"Lamb-Oseen: R={R:.0f} m,\n$\\Gamma$={G:.0f} m$^2$/s")
        ax[0].axvline(R, color=C[1], ls=":", lw=.8)
    ax[0].axhline(0, color="k", lw=.6)
    ax[0].set_xlabel("radius from fitted center (m)")
    ax[0].set_ylabel("velocity (m s$^{-1}$)")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Azimuthal velocity profile")
    ax[1].plot(rr, r["prof_zeta"], "o-", ms=3, lw=1, color=C[0], label="$\\zeta(r)$ local")
    if np.isfinite(r["oseen_radius"]):
        G, R = r["oseen_gamma"], r["oseen_radius"]
        ax[1].plot(rf, G / (np.pi * R ** 2) * np.exp(-(rf / R) ** 2), "-", lw=1.2,
                   color=C[1], label="Lamb-Oseen")
    ax[1].axhline(np.nanmedian(r["zeta_circ"]), color="0.4", ls="--", lw=.9,
                  label="cluster-mean $\\zeta$")
    ax[1].axhline(0, color="k", lw=.6)
    ax[1].set_xlabel("radius from fitted center (m)")
    ax[1].set_ylabel("$\\zeta$ (s$^{-1}$)")
    ax[1].legend(fontsize=7)
    ax[1].set_title("Local vorticity vs radius")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "radial_profile.png"))
    plt.close(fig)
    return ["vorticity.png", "kinematics.png", "tracks_center.png",
            "diagnostics.png", "radial_profile.png"]


def report(r):
    f = r["f"]
    z = r["zeta_circ"]
    ok = np.isfinite(z)
    med = np.nanmedian(z)
    mad = np.nanmedian(np.abs(z[ok] - med))
    lines = [
        f"window                : {_dt64(r['t'][0])}  ->  {_dt64(r['t'][-1])}",
        f"samples               : {len(r['t'])} cluster epochs, "
        f"{np.isfinite(r['u']).sum()} drifter-velocity values",
        f"Coriolis f at {r['lat0']:.3f}N  : {f:.3e} 1/s",
        "",
        "VORTICITY (circulation / Stokes, the primary estimator)",
        f"  median            : {med:+.3e} 1/s     Rossby {med/f:+.1f}",
        f"  MAD               : {mad:.3e} 1/s",
        f"  range (5-95%)     : {np.nanpercentile(z,5):+.3e} .. {np.nanpercentile(z,95):+.3e}",
        f"  rotation period   : {2*np.pi/abs(med/2)/3600:.2f} h  (as solid body, T = 4pi/|zeta|)",
        f"  leave-one-out 1sd : {np.nanmedian(r['zeta_jk_spread']):.3e} 1/s "
        f"({100*np.nanmedian(np.abs(r['zeta_jk_spread']/r['zeta_jk'])):.0f}% of |zeta|)",
        f"  LSQ formal 1sd    : {np.nanmedian(r['zeta_lsq_err']):.3e} 1/s "
        f"({100*np.nanmedian(r['zeta_lsq_err'])/abs(med):.0f}% of |zeta|)"
        "   <- optimistic, see notes",
        "",
        "CONSISTENCY",
        f"  constellation turned {r['turns']:+.2f} revolutions; "
        f"zeta/2 integrated predicts {np.nanmean(z)/2*(r['t'][-1]-r['t'][0])/1000/(2*np.pi):+.2f}",
        f"  divergence |delta|/|zeta| : contour {abs(np.nanmedian(r['delta_circ']))/abs(med):.3f}, "
        f"area {abs(np.nanmedian(r['delta_area']))/abs(med):.3f}, "
        f"LSQ {abs(np.nanmedian(r['delta_lsq']))/abs(med):.3f}",
        f"  Okubo-Weiss < 0 in  {100*np.nanmean(r['OW']<0):.0f}% of windows "
        "(rotation dominates strain)",
        f"  |strain|/|zeta|     : {np.nanmedian(np.hypot(r['sig_n'],r['sig_s']))/abs(med):.2f}",
        "",
        "CLUSTER GEOMETRY",
        f"  area   median {np.nanmedian(r['area'])/1e6:.3f} km^2   "
        f"({np.nanmin(r['area'])/1e6:.3f} .. {np.nanmax(r['area'])/1e6:.3f})",
        f"  scale  sqrt(area) median {np.sqrt(np.nanmedian(r['area'])):.0f} m",
        f"  aspect median {np.nanmedian(r['aspect']):.2f}, "
        f"p95 {np.nanpercentile(r['aspect'],95):.2f}",
        f"  eddy radius 2|v|/|zeta| ~ "
        f"{2*np.nanmedian(np.hypot(r['u'],r['v']))/abs(med):.0f} m",
        f"  fit condition number median {np.nanmedian(r['cond']):.2f}, "
        f"residual {np.nanmedian(r['rms']):.4f} m/s",
        "",
        "VORTEX STRUCTURE (Lamb-Oseen fit to the radial profile)",
        f"  core radius R      : {r['oseen_radius']:.0f} m",
        f"  circulation Gamma  : {r['oseen_gamma']:.0f} m^2/s",
        f"  core vorticity     : {r['oseen_zeta_core']:+.3e} 1/s  "
        f"(Rossby {r['oseen_zeta_core']/f:+.0f})",
        f"  |v_r|/|v_theta|    : {np.nanmedian(np.abs(r['prof_vr']))/np.nanmedian(np.abs(r['prof_vt'])):.3f}"
        "   (small => coherent, not dispersing)",
    ]
    return "\n".join(lines)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("-i", "--input", default=os.path.join(DATA, "drifters"),
                   help="directory holding mwb*_gps_timeseries.nc")
    p.add_argument("-o", "--output", default="eddy_out",
                   help="output directory (regenerable; not in git)")
    p.add_argument("--start", default="2023-05-22T05:10:54")
    p.add_argument("--end", default="2023-05-23T06:34:00")
    p.add_argument("--blocks", type=int, default=4,
                   help="block averages per 2048-sample burst (default 4 = 256 s)")
    p.add_argument("--window", type=float, default=1800.0,
                   help="least-squares window, seconds")
    p.add_argument("--step", type=float, default=600.0)
    p.add_argument("--no-figures", action="store_true")
    a = p.parse_args(argv)

    t0 = np.datetime64(a.start).astype("datetime64[ms]").astype(np.int64)
    t1 = np.datetime64(a.end).astype("datetime64[ms]").astype(np.int64)
    r = run(a.input, a.output, t0, t1, blocks=a.blocks,
            window_s=a.window, step_s=a.step)
    print(report(r))
    if not a.no_figures:
        for f in figures(r, a.output):
            print(f"  wrote {os.path.join(a.output, f)}")
    print(f"  wrote {os.path.join(a.output, 'eddy_kinematics.npz')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
