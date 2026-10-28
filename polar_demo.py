"""
One-off: show bn180715741 in the current flat projection next to a polar one.

The burst sits at Dec +79.8 and both 90% regions are caps that contain the
north celestial pole. A region containing a pole has no closed boundary in
RA/Dec, so the flat panel draws only the cap's lower rim and the cross looks
stranded outside it. The polar panel puts the pole at the centre, where the
same contour closes into the circle it actually is.

Not wired into the viewer -- this is a look at whether the reprojection is
worth doing.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TRIGGER = "bn180715741"
SRC = f"site/vec/bursts/{TRIGGER}.json"
OUT = f"/tmp/claude-1002/-home-aakyuz-Fermi/47126f7a-c737-4295-ba11-7d3148306212/scratchpad/polar_{TRIGGER}.png"

OFFICIAL_C = "#2f6fd0"
T90_C = "#b8690a"
PEAK_C = "#d9534f"


def paths(rec, key):
    for flat in rec.get(key, []):
        yield np.array(flat[0::2], dtype=float), np.array(flat[1::2], dtype=float)


def main():
    with open(SRC) as handle:
        rec = json.load(handle)

    layers = [
        ("O", OFFICIAL_C, "official Fermi GBM", rec["oc"], "-", 2.0),
        ("A", T90_C, "mine, full T90", rec["ac"], "--", 2.0),
        ("P", PEAK_C, "mine, best S/N", rec["pc"], "-", 1.6),
    ]

    fig = plt.figure(figsize=(12.4, 5.9), dpi=110)

    # ---------------------------------------------------- flat RA/Dec (now)
    axf = fig.add_subplot(1, 2, 1)
    centre = rec["oc"][0]
    for key, colour, label, cen, style, width in layers:
        for i, (ra, dec) in enumerate(paths(rec, key)):
            rw = centre + (ra - centre + 180.0) % 360.0 - 180.0
            axf.plot(rw, dec, color=colour, lw=width, ls=style,
                     label=label if i == 0 else None)
        rw = centre + (cen[0] - centre + 180.0) % 360.0 - 180.0
        axf.plot([rw], [cen[1]], marker="P" if key == "O" else "X",
                 color=colour, ms=11, mec="white", mew=1.1, ls="none")
    axf.invert_xaxis()
    axf.set_xlabel("RA (deg)")
    axf.set_ylabel("Dec (deg)")
    axf.grid(alpha=.25, lw=.6)
    axf.set_title("Flat RA/Dec — what the viewer draws now\n"
                  "the official rim sits far below its own cross", fontsize=10)
    axf.legend(fontsize=8, loc="lower left", framealpha=.92)

    # ------------------------------------------------ polar, pole at centre
    axp = fig.add_subplot(1, 2, 2, projection="polar")
    axp.set_theta_zero_location("S")
    axp.set_theta_direction(1)          # RA increases counterclockwise
    for key, colour, label, cen, style, width in layers:
        for i, (ra, dec) in enumerate(paths(rec, key)):
            th = np.radians(ra % 360.0)
            r = 90.0 - dec
            # close the loop: a cap's boundary wraps all the way round in RA
            th = np.append(th, th[0])
            r = np.append(r, r[0])
            axp.plot(th, r, color=colour, lw=width, ls=style,
                     label=label if i == 0 else None)
            axp.fill(th, r, color=colour, alpha=.10, lw=0)
        axp.plot([np.radians(cen[0] % 360.0)], [90.0 - cen[1]],
                 marker="P" if key == "O" else "X", color=colour, ms=11,
                 mec="white", mew=1.1, ls="none")

    axp.set_rmin(0)
    axp.set_rmax(54)
    axp.set_rgrids([10, 20, 30, 40, 50],
                   labels=["+80", "+70", "+60", "+50", "+40"], fontsize=7.5,
                   angle=112)
    axp.set_thetagrids(range(0, 360, 45),
                       labels=[f"{d}h" for d in range(0, 24, 3)], fontsize=7.5)
    axp.grid(alpha=.3, lw=.6)
    axp.set_title("Polar, north pole at centre — the fix\n"
                  "the same contours close into caps around the cross",
                  fontsize=10)
    axp.legend(fontsize=8, loc="upper left", bbox_to_anchor=(-0.16, 1.06),
               framealpha=.92)
    # The ring labels already read as declinations; spelling that out beside
    # them collided with the 9h spoke, so it goes under the axes instead.
    axp.text(0.5, -0.09, "rings are declination; centre of the plot = pole, +90",
             transform=axp.transAxes, fontsize=7.5, color="#6b7686", ha="center")

    fig.suptitle(f"{TRIGGER} — T90 {rec['t90']} s, centroid Dec "
                 f"{rec['oc'][1]:+.1f}, official 90% area {rec['oa']} sq deg",
                 fontsize=11)
    plt.tight_layout(rect=(0, 0, 1, 0.94))
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    plt.savefig(OUT, bbox_inches="tight")
    plt.close("all")
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1024:.0f} KB)")


if __name__ == "__main__":
    main()
