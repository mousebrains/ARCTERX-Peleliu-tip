#! /usr/bin/env python3

"""One-page summary figure of the Peleliu tip vortex analysis.

Reads eddy_out/eddy_kinematics.npz (written by eddy_analysis.py) and renders a
presentation-quality panel: headline numbers, the vorticity time series with its
leave-one-out uncertainty, the radial vorticity structure, and the drifter
tracks in the eddy frame.

    python3 eddy_summary_figure.py [--input eddy_out] [--dark]
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

# Validated categorical palette (see dataviz reference palette; the three-slot
# subset clears the all-pairs CVD and normal-vision floors in both modes).
LIGHT = dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", ink3="#8a8880",
             grid="#e3e2dd", s1="#2a78d6", s2="#eb6834", s3="#1baf7a")
DARK = dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", ink3="#8a8880",
            grid="#333331", s1="#3987e5", s2="#d95926", s3="#199e70")


def _tile(ax, value, label, sub, P, color=None):
    ax.axis("off")
    ax.text(0, .78, value, fontsize=19, fontweight="600", va="top",
            color=color or P["ink"], transform=ax.transAxes)
    ax.text(0, .34, label, fontsize=8.5, va="top", color=P["ink2"],
            transform=ax.transAxes)
    ax.text(0, .10, sub, fontsize=7.5, va="top", color=P["ink3"],
            transform=ax.transAxes)


def _clean(ax, P):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(P["grid"])
    ax.tick_params(colors=P["ink2"], labelsize=7.5, length=3, width=.8)
    ax.grid(True, color=P["grid"], lw=.6, alpha=.9)
    ax.set_axisbelow(True)
    for lb in list(ax.get_xticklabels()) + list(ax.get_yticklabels()):
        lb.set_color(P["ink2"])


def build(npz, out, dark=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    P = DARK if dark else LIGHT
    r = np.load(npz, allow_pickle=True)
    f = float(r["f"])
    T = r["t"].astype(np.int64).astype("datetime64[ms]").astype("datetime64[s]").tolist()
    z, sp = r["zeta_circ"], r["zeta_jk_spread"]
    zj = r["zeta_jk"]
    med = np.nanmedian(z)

    fig = plt.figure(figsize=(11.5, 8.2), facecolor=P["surface"])
    gs = fig.add_gridspec(3, 4, height_ratios=[.42, 1.15, 1.0],
                          hspace=.62, wspace=.34,
                          left=.075, right=.928, top=.895, bottom=.115)

    fig.text(.065, .965, "Peleliu tip vortex", fontsize=16, fontweight="600",
             color=P["ink"])
    fig.text(.065, .933,
             "Relative vorticity from four CORDC mini wave buoys  ·  "
             "2023-05-22 05:11 to 05-23 06:34 UTC  ·  25.4 h",
             fontsize=8.5, color=P["ink2"])

    # ---- headline numbers -------------------------------------------------
    _tile(fig.add_subplot(gs[0, 0]), f"{med:.2e}".replace("e-0", "×10⁻") + " s⁻¹",
          "median relative vorticity ζ", "circulation / Stokes, ±17 %", P, P["s1"])
    _tile(fig.add_subplot(gs[0, 1]), f"{med/f:.0f}", "Rossby number  ζ/f",
          f"f = {f:.2e} s⁻¹ at 6.99°N", P)
    _tile(fig.add_subplot(gs[0, 2]), f"{2*np.pi/abs(med/2)/3600:.1f} h",
          "core rotation period", "vortex stayed coherent throughout", P)
    _tile(fig.add_subplot(gs[0, 3]), f"{float(r['oseen_radius']):.0f} m",
          "Lamb–Oseen core radius", f"Γ = {float(r['oseen_gamma']):.0f} m² s⁻¹", P)

    # ---- vorticity time series -------------------------------------------
    ax = fig.add_subplot(gs[1, :])
    _clean(ax, P)
    ax.fill_between(T, (zj - sp) * 1e3, (zj + sp) * 1e3, color=P["s1"],
                    alpha=.16, lw=0, zorder=1)
    ax.plot(T, 2e3 * r["omega_pos"], color=P["s3"], lw=1.0, alpha=.85, zorder=2)
    FT = r["fit_t"].astype(np.int64).astype("datetime64[ms]").astype("datetime64[s]").tolist()
    ax.plot(FT, 1e3 * r["zeta_lsq"], color=P["s2"], lw=1.4, zorder=3)
    ax.plot(T, 1e3 * z, color=P["s1"], lw=1.8, zorder=4)
    ax.axhline(0, color=P["ink3"], lw=.8, zorder=1)
    ax.axhline(1e3 * med, color=P["s1"], lw=.8, ls=(0, (4, 3)), alpha=.55, zorder=1)
    ax.set_ylim(-4.2, 1.4)
    ax.set_ylabel("ζ   (10$^{-3}$ s$^{-1}$)", fontsize=9, color=P["ink2"],
                  labelpad=6)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d %b\n%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    sec = ax.secondary_yaxis("right", functions=(lambda a: a * 1e-3 / f,
                                                 lambda a: a * f * 1e3))
    sec.set_ylabel("Rossby number  ζ/f", fontsize=9, color=P["ink2"],
                   labelpad=8)
    sec.tick_params(colors=P["ink2"], labelsize=7.5, length=3, width=.8)
    sec.spines["right"].set_color(P["grid"])
    # Legend sits in its own strip above the axes: inline labels collided with
    # the traces.  Swatch + text keeps identity off color alone, which also
    # satisfies the relief rule for the aqua step on the light surface.
    ax.set_title("Vorticity is anticyclonic and sustained for the full record; "
                 "shading is the leave-one-drifter-out spread",
                 fontsize=9.5, color=P["ink2"], loc="left", pad=26)
    for xf, txt, col in ((.000, "circulation / Stokes  (primary)", P["s1"]),
                         (.285, "least-squares gradient", P["s2"]),
                         (.505, "constellation rotation (positions only)", P["s3"])):
        ax.annotate("", xy=(xf + .026, 1.035), xytext=(xf, 1.035),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="-", color=col, lw=2.4))
        ax.text(xf + .034, 1.035, txt, transform=ax.transAxes, fontsize=8,
                color=P["ink2"], va="center")

    # ---- radial structure -------------------------------------------------
    ax = fig.add_subplot(gs[2, :2])
    _clean(ax, P)
    rr, pz = r["prof_r"], r["prof_zeta"] * 1e3
    good = np.isfinite(pz) & (r["prof_n"] > 30)
    G, R = float(r["oseen_gamma"]), float(r["oseen_radius"])
    rf = np.linspace(rr[good].min(), rr[good].max(), 200)
    ax.plot(rf, 1e3 * G / (np.pi * R ** 2) * np.exp(-(rf / R) ** 2),
            color=P["s2"], lw=1.6, zorder=2)
    ax.plot(rr[good], pz[good], "o-", color=P["s1"], lw=1.6, ms=6,
            mec=P["surface"], mew=1.4, zorder=3)
    ax.axhline(1e3 * med, color=P["ink3"], lw=.9, ls=(0, (4, 3)), zorder=1)
    ax.annotate("cluster-mean ζ", xy=(rr[good].max(), 1e3 * med), fontsize=7.5,
                color=P["ink3"], va="bottom", ha="right")
    ax.annotate("observed ζ(r)", xy=(rr[good][0], pz[good][0]),
                xytext=(26, 4), textcoords="offset points",
                fontsize=8, fontweight="600", color=P["s1"])
    ax.annotate("Lamb–Oseen fit", xy=(rf[40], 1e3 * G / (np.pi * R ** 2)
                                      * np.exp(-(rf[40] / R) ** 2)),
                xytext=(10, -20), textcoords="offset points",
                fontsize=8, fontweight="600", color=P["s2"],
                arrowprops=dict(arrowstyle="-", color=P["s2"], lw=.8,
                                shrinkA=2, shrinkB=2))
    ax.set_xlabel("radius from fitted eddy center  (m)", fontsize=9, color=P["ink2"])
    ax.set_ylabel("ζ   (10$^{-3}$ s$^{-1}$)", fontsize=9, color=P["ink2"],
                  labelpad=6)
    ax.set_title("The vortex is not solid-body: ζ falls off with radius",
                 fontsize=9.5, color=P["ink2"], loc="left", pad=8)

    # ---- eddy-frame tracks ------------------------------------------------
    ax = fig.add_subplot(gs[2, 2:])
    _clean(ax, P)
    ok = np.isfinite(r["center_x"])
    cxi = np.interp(r["t"], r["fit_t"][ok], r["center_x"][ok])
    cyi = np.interp(r["t"], r["fit_t"][ok], r["center_y"][ok])
    for j in range(r["x"].shape[1]):
        ax.plot((r["x"][:, j] - cxi) / 1e3, (r["y"][:, j] - cyi) / 1e3,
                lw=.75, color=P["s1"], alpha=.5)
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(R / 1e3 * np.cos(th), R / 1e3 * np.sin(th), color=P["s2"], lw=1.6)
    ax.annotate(f"core R = {R:.0f} m", xy=(R / 1e3 * .7, R / 1e3 * .72),
                fontsize=8, fontweight="600", color=P["s2"])
    ax.plot(0, 0, "+", ms=11, color=P["ink"], mew=1.6)
    ax.set_aspect("equal")
    ax.set_xlabel("east of center  (km)", fontsize=9, color=P["ink2"])
    ax.set_ylabel("north of center  (km)", fontsize=9, color=P["ink2"])
    ax.set_title("Drifter tracks in the eddy frame (all four)",
                 fontsize=9.5, color=P["ink2"], loc="left", pad=8)

    fig.text(.075, .062,
             "Vorticity from circulation around the drifter polygon (Stokes), cross-checked against a "
             "least-squares velocity-gradient fit and a positions-only rotation rate.\n"
             "Divergence is < 2 % of |ζ| by three independent estimators; Okubo–Weiss < 0 in 100 % of "
             "windows.  Each ζ is an area average over the polygon.\n"
             "Shaded band and the ±17 % are leave-one-drifter-out spreads — field curvature, not noise: "
             "they do not shrink with more averaging.",
             fontsize=7, color=P["ink3"], va="top", linespacing=1.7)

    for a in fig.get_axes():
        a.set_facecolor(P["surface"])
    fig.savefig(out, dpi=200, facecolor=P["surface"])
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("-i", "--input", default="eddy_out",
                   help="directory holding eddy_kinematics.npz")
    p.add_argument("--dark", action="store_true")
    a = p.parse_args(argv)
    npz = os.path.join(a.input, "eddy_kinematics.npz")
    out = os.path.join(a.input, "summary_dark.png" if a.dark else "summary.png")
    print("wrote", build(npz, out, dark=a.dark))
    return 0


if __name__ == "__main__":
    sys.exit(main())
