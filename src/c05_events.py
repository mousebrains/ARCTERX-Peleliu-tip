"""C05 rotation, principal axis, subtidal regimes and the high-flow event list.

Regenerates data/adcp/c05_high_flow_events.csv from data/adcp/c05_depth_avg.npz
(stored WITHOUT the compass rotation) for a chosen absolute rotation.

    python3 src/c05_events.py                      # the adopted rotation (19.5 deg)
    python3 src/c05_events.py --rotation 14.4      # the superseded one, for comparison
    python3 src/c05_events.py --write              # also rewrite the CSV

The rotation is W_corr = W * exp(-i * rotation).  C05 reads 21.09 deg counterclockwise of
the Sig1000 (complex transfer, WAMOS; reproduced here as 21.00 deg).  The Sig1000 in turn
reads 1.6 deg clockwise of the R/V Thompson wh300 over the crest (160 pairs within 0.3 km,
Hydrographer-Analysis notes/10 section 2; MAD 2.1 deg), so C05 reads 21.09 - 1.6 = 19.5
deg counterclockwise of true.  The earlier 14.4 deg used a 6.6 deg Signature offset from a
network fit that predates the Thompson cross-check (Pat Welch, 2026-10-05).

Definitions (match PRESSURE_ANALYSIS.md section 4):
  principal axis   eigenvector of the east/north covariance, as a bearing in [90, 270)
                   for the positive end; ellipticity = sqrt(minor/major eigenvalue)
  along            velocity on that axis, + toward the axis bearing
  subtidal         daily means (UTC days) of the along-axis velocity
  events           contiguous runs of |along| above its 85th percentile lasting > 1 h;
                   direction from the sign of the along-axis velocity at the peak
"""
import argparse
import datetime as dt
import os

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ADOPTED = 19.5


def load():
    d = np.load(f"{HERE}/data/adcp/c05_depth_avg.npz")
    return d["time_ms"] / 1000.0, d["east"].astype(float) + 1j * d["north"].astype(float)


def axis(w):
    C = np.cov(np.vstack([w.real, w.imag]))
    ev, V = np.linalg.eigh(C)
    major = V[:, 1]
    bearing = np.degrees(np.arctan2(major[0], major[1])) % 360
    if not 90 <= bearing < 270:
        bearing = (bearing + 180) % 360
    return bearing, np.sqrt(ev[0] / ev[1])


def analyse(rotation, write=False):
    t, w = load()
    wc = w * np.exp(-1j * np.radians(rotation))
    b, ell = axis(wc)
    u = np.array([np.sin(np.radians(b)), np.cos(np.radians(b))])
    along = wc.real * u[0] + wc.imag * u[1]
    day = np.floor(t / 86400).astype(int)
    dmean = {k: along[day == k].mean() for k in np.unique(day)}
    sub = np.array([dmean[k] for k in day])
    agree = np.mean(np.sign(along) == np.sign(sub))
    thr = np.percentile(np.abs(along), 85)
    hi = np.abs(along) > thr
    events = []
    i = 0
    while i < len(hi):
        if hi[i]:
            j = i
            while j + 1 < len(hi) and hi[j + 1] and t[j + 1] - t[j] < 1800:
                j += 1
            if t[j] - t[i] > 3600:
                k = i + np.argmax(np.abs(along[i:j + 1]))
                events.append((t[i], t[j], (t[j] - t[i]) / 3600, along[k] * 100, "EAST" if along[k] > 0 else "WEST"))
            i = j + 1
        else:
            i += 1
    ne = sum(e[4] == "EAST" for e in events)
    print(f"rotation {rotation:.1f} deg: principal axis {b:.1f}/{(b + 180) % 360:.1f} deg true, ellipticity {ell:.2f}")
    print(f"  sd east {wc.real.std()*100:.1f}, north {wc.imag.std()*100:.1f} cm/s; threshold {thr*100:.1f} cm/s")
    print(f"  sign(total) = sign(subtidal daily mean) {agree*100:.0f} %")
    print(f"  events > 1 h: {len(events)} ({ne} east, {len(events) - ne} west)")
    may = [e for e in events if dt.datetime.fromtimestamp(e[0], dt.timezone.utc).strftime("%Y-%m") == "2023-05"]
    for e in may:
        if dt.datetime.fromtimestamp(e[0], dt.timezone.utc).day in (21, 22, 23):
            print(f"    {dt.datetime.fromtimestamp(e[0], dt.timezone.utc):%m-%d %H:%M} - "
                  f"{dt.datetime.fromtimestamp(e[1], dt.timezone.utc):%m-%d %H:%M}  peak {e[3]:+.0f} cm/s {e[4]}")
    print("  May 2023 daily subtidal (cm/s):", " ".join(
        f"{dt.datetime.fromtimestamp(k*86400, dt.timezone.utc):%d}:{dmean[k]*100:+.0f}" for k in sorted(dmean)
        if dt.datetime.fromtimestamp(k*86400, dt.timezone.utc).strftime("%Y-%m") == "2023-05"
        and dt.datetime.fromtimestamp(k*86400, dt.timezone.utc).day >= 14))
    if write:
        fn = f"{HERE}/data/adcp/c05_high_flow_events.csv"
        with open(fn, "w") as fp:
            fp.write(f'"# C05 depth-averaged current, absolute rotation exp(-i*{rotation:.1f}deg) applied (src/c05_events.py)"\n')
            fp.write(f'"# along principal axis {b:.1f}/{(b + 180) % 360:.1f} deg true; +=toward {b:.0f}, -=toward {(b + 180) % 360:.0f}"\n')
            fp.write("start_utc,end_utc,duration_h,peak_along_cm_s,direction\n")
            for e in events:
                fp.write(f"{dt.datetime.fromtimestamp(e[0], dt.timezone.utc):%Y-%m-%dT%H:%M:%S},"
                         f"{dt.datetime.fromtimestamp(e[1], dt.timezone.utc):%Y-%m-%dT%H:%M:%S},"
                         f"{e[2]:.2f},{e[3]:.1f},{e[4]}\n")
        print("  wrote", fn)
    return dict(axis=b, ell=ell, along=along, events=events)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rotation", type=float, default=ADOPTED, help="absolute rotation, deg (default %(default)s)")
    ap.add_argument("--write", action="store_true", help="rewrite data/adcp/c05_high_flow_events.csv")
    a = ap.parse_args()
    r = analyse(a.rotation, a.write)
    if a.rotation != 14.4:
        old = analyse(14.4)
        flip = np.mean(np.sign(old["along"]) != np.sign(r["along"]))
        print(f"  along-axis sign changes vs 14.4 deg: {flip*100:.2f} % of ensembles")


if __name__ == "__main__":
    main()
