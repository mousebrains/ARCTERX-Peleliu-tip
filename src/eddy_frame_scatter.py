#! /usr/bin/env python3

"""Eddy-frame position of each drifter, colored by hours from the common start.

One panel per drifter (small multiples) sharing axes, scale and color ramp, so
the four are directly comparable: any difference between panels is a difference
in the drifters, not in the plotting.

    python3 eddy_frame_scatter.py [--input eddy_out] [--dark]

Color encodes elapsed time, which is a magnitude, so the ramp is sequential --
a single hue running light to dark.  A rainbow (jet/turbo) would imply
categories where there are none and is not perceptually uniform.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

LIGHT = dict(surface="#fcfcfb", ink="#0b0b0b", ink2="#52514e", ink3="#8a8880",
             grid="#e3e2dd", ramp=("#d3e3f7", "#0f3a68"))
DARK = dict(surface="#1a1a19", ink="#ffffff", ink2="#c3c2b7", ink3="#8a8880",
            grid="#333331", ramp=("#17324f", "#8fc2f5"))

DRIFTERS = ("mwb458d02", "mwb788d01", "mwb790d01", "mwb793d02")


def build(npz, out, dark=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize

    P = DARK if dark else LIGHT
    r = np.load(npz, allow_pickle=True)
    cmap = LinearSegmentedColormap.from_list("elapsed", list(P["ramp"]), N=256)

    # eddy-frame coordinates: drifter position minus the fitted center
    ok = np.isfinite(r["center_x"])
    cx = np.interp(r["t"], r["fit_t"][ok], r["center_x"][ok])
    cy = np.interp(r["t"], r["fit_t"][ok], r["center_y"][ok])
    X = (r["x"] - cx[:, None]) / 1e3
    Y = (r["y"] - cy[:, None]) / 1e3
    hours = (r["t"] - r["t"][0]) / 3.6e6
    norm = Normalize(0, hours.max())

    # generous enough that nothing is silently clipped; the count of any
    # points outside is asserted below rather than hidden
    lim = np.nanpercentile(np.abs(np.concatenate([X.ravel(), Y.ravel()])), 99.0)
    lim = np.ceil(lim * 4) / 4
    nout = int(np.nansum((np.abs(X) > lim) | (np.abs(Y) > lim)))

    fig, axes = plt.subplots(2, 2, figsize=(9.4, 9.0), facecolor=P["surface"],
                             sharex=True, sharey=True)
    fig.subplots_adjust(left=.085, right=.87, top=.885, bottom=.145,
                        wspace=.12, hspace=.16)

    fig.text(.085, .955, "Drifter position in the eddy reference frame",
             fontsize=15, fontweight="600", color=P["ink"])
    fig.text(.085, .925,
             "Each panel is one drifter, relative to the fitted vortex center  ·  "
             "color is hours from the common start (2023-05-22 05:14 UTC)",
             fontsize=8.5, color=P["ink2"])

    R = float(r["oseen_radius"]) / 1e3
    th = np.linspace(0, 2 * np.pi, 240)
    for k, ax in enumerate(axes.ravel()):
        ax.set_facecolor(P["surface"])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(P["grid"])
        ax.grid(True, color=P["grid"], lw=.6)
        ax.set_axisbelow(True)
        ax.tick_params(colors=P["ink2"], labelsize=8, length=3, width=.8)

        # core circle and center marker give every panel the same reference
        ax.plot(R * np.cos(th), R * np.sin(th), color=P["ink3"], lw=1.0,
                ls=(0, (5, 4)), zorder=2)
        ax.plot(0, 0, "+", ms=10, color=P["ink2"], mew=1.4, zorder=3)

        g = np.isfinite(X[:, k]) & np.isfinite(Y[:, k])
        sc = ax.scatter(X[g, k], Y[g, k], c=hours[g], cmap=cmap, norm=norm,
                        s=11, lw=0, zorder=4)
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_title(DRIFTERS[k], fontsize=10, fontweight="600",
                     color=P["ink"], loc="left", pad=6)
        if k == 0:
            ax.annotate(f"core R = {R*1e3:.0f} m",
                        xy=(R * np.cos(np.pi / 4), R * np.sin(np.pi / 4)),
                        xytext=(10, 8), textcoords="offset points",
                        fontsize=7.5, color=P["ink3"])
        if k >= 2:
            ax.set_xlabel("east of eddy center  (km)", fontsize=9, color=P["ink2"])
        if k % 2 == 0:
            ax.set_ylabel("north of eddy center  (km)", fontsize=9, color=P["ink2"])

    cax = fig.add_axes([.895, .145, .018, .74])
    cb = fig.colorbar(sc, cax=cax)
    cb.set_label("hours from common start", fontsize=9, color=P["ink2"])
    cb.ax.tick_params(colors=P["ink2"], labelsize=8, length=3, width=.8)
    cb.outline.set_edgecolor(P["grid"])

    fig.text(.085, .082,
             "Dashed circle is the Lamb–Oseen core radius fitted to the pooled radial "
             "profile.  Center is the elliptic critical point of the windowed\n"
             "velocity-gradient fit and is conditional on the eddy translation "
             "estimate, so a slow drift of the cloud off the origin is a center "
             f"error,\nnot a drifter leaving the vortex.  {nout} of {X.size} points "
             "fall outside the axes.",
             fontsize=7, color=P["ink3"], va="top", linespacing=1.7)

    fig.savefig(out, dpi=200, facecolor=P["surface"])
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("-i", "--input", default="eddy_out",
                   help="directory holding eddy_kinematics.npz")
    p.add_argument("--dark", action="store_true")
    a = p.parse_args(argv)
    npz = os.path.join(a.input, "eddy_kinematics.npz")
    out = os.path.join(a.input,
                       "eddy_frame_scatter_dark.png" if a.dark
                       else "eddy_frame_scatter.png")
    print("wrote", build(npz, out, dark=a.dark))
    return 0


if __name__ == "__main__":
    sys.exit(main())
