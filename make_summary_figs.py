"""
Summary figures pooling both samples: the window change, and the energy sweep.
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SITEDIR = "site"
FIGDIR = os.path.join(SITEDIR, "figures")
MINE = "#d9534f"
OFFICIAL = "#2f6fd0"
GOOD = "#2e7d55"
CRIT = "#bf3b36"


def load(name):
    with open(os.path.join(SITEDIR, name)) as handle:
        return json.load(handle)


def paired():
    """(trigger, t90_sep, peak_sep, grb_type, sample) for both samples."""
    out = []
    a1 = {r["trigger"]: r for r in load("results.json")["results"]}
    b1 = {r["trigger"]: r for r in load("results_peak.json")["results"]}
    for k in sorted(set(a1) & set(b1)):
        out.append((k, a1[k]["sep_mine_official"], b1[k]["sep_mine_official"],
                    b1[k]["grb_type"], "A"))
    a2 = {r["trigger"]: r for r in load("holdout_t90.json")["results"]}
    b2 = {r["trigger"]: r for r in load("holdout_50_300.json")["results"]}
    for k in sorted(set(a2) & set(b2)):
        out.append((k, a2[k]["sep_mine_official"], b2[k]["sep_mine_official"],
                    b2[k]["grb_type"], "B"))
    return out


def slope_chart(rows):
    """One line per burst, T90 window on the left, best-S/N on the right."""
    fig, ax = plt.subplots(figsize=(7.4, 5.4), dpi=110)

    for _, before, after, gtype, _sample in rows:
        better = after < before - 1e-9
        ax.plot([0, 1], [before, after],
                color=(GOOD if better else CRIT),
                alpha=.75, lw=1.4, zorder=2,
                marker="o" if gtype == "long" else "s", ms=4.5,
                markerfacecolor="white", markeredgewidth=1.3)

    before = np.array([r[1] for r in rows])
    after = np.array([r[2] for r in rows])
    ax.plot([0, 1], [np.median(before), np.median(after)],
            color="#141b26", lw=3.4, zorder=4, marker="D", ms=8,
            label=f"median  {np.median(before):.2f}$\\degree$ "
                  f"$\\rightarrow$ {np.median(after):.2f}$\\degree$")

    ax.axhspan(0, 4.2, color=OFFICIAL, alpha=.07, zorder=0)
    ax.text(1.32, 2.0, "GBM systematic\ncore (2.3-4.2$\\degree$)",
            fontsize=7.5, color=OFFICIAL, va="center")

    ax.set_yscale("log")
    ax.set_xlim(-0.33, 1.92)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["full T90 window", "best-S/N window"], fontsize=10)
    ax.set_ylabel("separation from official map (deg)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", alpha=.25, lw=.6)
    ax.set_title(f"Effect of the source window, {len(rows)} bursts "
                 f"(circles long, squares short)", fontsize=9.5)
    ax.legend(fontsize=8, loc="upper right", framealpha=.92)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "summary_window.png"), bbox_inches="tight")
    plt.close("all")


def band_chart():
    bands = load("sweep_main_summary.json")["bands"]
    bands = sorted(bands, key=lambda b: b["median"])
    labels = [b["label"].split("(")[0].strip() for b in bands]
    medians = [b["median"] for b in bands]
    is_std = [(b["lo"], b["hi"]) == (50, 300) for b in bands]

    fig, ax = plt.subplots(figsize=(7.4, 3.9), dpi=110)
    ypos = np.arange(len(bands))
    ax.barh(ypos, medians,
            color=[OFFICIAL if std else "#b9c2d0" for std in is_std],
            height=.66, zorder=3)
    for y, value, std in zip(ypos, medians, is_std):
        ax.text(value + .12, y, f"{value:.2f}" + ("  standard" if std else ""),
                va="center", fontsize=8,
                color=OFFICIAL if std else "#46525f",
                fontweight="bold" if std else "normal")

    ax.set_yticks(ypos)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("median separation from official map (deg)", fontsize=9)
    ax.set_xlim(0, max(medians) * 1.26)
    ax.tick_params(labelsize=8)
    ax.grid(axis="x", alpha=.25, lw=.6, zorder=0)
    ax.set_title("Localization band: narrowing costs accuracy, widening gains nothing",
                 fontsize=9.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "summary_bands.png"), bbox_inches="tight")
    plt.close("all")


if __name__ == "__main__":
    os.makedirs(FIGDIR, exist_ok=True)
    rows = paired()
    slope_chart(rows)
    band_chart()
    print(f"wrote summary_window.png ({len(rows)} bursts) and summary_bands.png")
