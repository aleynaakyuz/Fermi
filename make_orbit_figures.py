"""
Figures explaining Fermi's orbit and how it meets the South Atlantic Anomaly.

Everything else on the SAA page is hand-written SVG, because it is all line
charts over data this project computed. These are maps, and maps want real
coastlines, so they are rendered with cartopy and embedded as PNG.

Three figures:
    orbit_explained   the zone on a world map, five consecutive orbits
                      stepping west across it, and one day as a timeline
    orbit_all_zoom    every orbit of one day over the anomaly, each crossing
                      drawn in its own orbit's colour
    orbit_all_small   the same day as small multiples, one orbit per panel

Two choices worth knowing about. Orbits are split at the DESCENDING node,
which falls on the far side of the Earth from the anomaly: splitting at the
ascending node cuts a crossing in half and makes the per-orbit durations
disagree with the timeline. And each crossing is coloured by its orbit rather
than a single red -- ten overlapping red segments are an unreadable
crosshatch, which is the whole problem these figures exist to avoid.

Usage:
    python make_orbit_figures.py --day 220827
"""
import argparse
import glob
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np

DAILY = "data/daily"
OUTDIR = "site/figures"

INK = "#141b26"
MUT = "#6b7686"
HOT = "#d9534f"
RULE = "#d8dee8"


def load(day):
    from astropy.io import fits
    path = glob.glob(os.path.join(DAILY, day, "glg_poshist_all_*.fit"))[0]
    with fits.open(path) as handle:
        data = handle[1].data
        lon = np.asarray(data["SC_LON"], dtype=float)
        return (np.asarray(data["SCLK_UTC"], dtype=float),
                np.asarray(data["SC_LAT"], dtype=float),
                np.where(lon > 180, lon - 360, lon),
                ((np.asarray(data["FLAGS"]).astype(np.uint16) >> 1) & 1).astype(bool))


def saa_hull():
    """Convex hull of every flagged position across every day on disk."""
    from shapely.geometry import MultiPoint
    lats, lons = [], []
    for day in sorted(os.listdir(DAILY)):
        try:
            _, la, lo, sa = load(day)
        except Exception:
            continue
        lats.append(la[sa])
        lons.append(lo[sa])
    pts = np.column_stack([np.concatenate(lons), np.concatenate(lats)])
    return np.array(MultiPoint(pts).convex_hull.exterior.coords)


def split_seam(lon):
    """Indices grouped so no drawn segment jumps the antimeridian.

    Longitude is a wrapped coordinate: two samples a second apart can read
    +179.96 then -179.98, eight kilometres of motion and a 360 degree jump in
    the number. Joined literally that draws a horizontal line racing back
    across the whole map.
    """
    return np.split(np.arange(len(lon)),
                    np.flatnonzero(np.abs(np.diff(lon)) > 180) + 1)


def orbit_slices(lat):
    """One slice per orbit, cut at the descending node."""
    node = np.flatnonzero((lat[:-1] > 0) & (lat[1:] <= 0))
    return [slice(node[k], node[k + 1]) for k in range(len(node) - 1)]


def basemap(ax, extent, label_size=9):
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.OCEAN, facecolor="#f7fafd")
    ax.add_feature(cfeature.LAND, facecolor="#e9edf3",
                   edgecolor="#b9c3d0", linewidth=.7)
    gl = ax.gridlines(draw_labels=True, linewidth=.4, color="#e3e8ef")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": label_size, "color": MUT}
    gl.ylabel_style = {"size": label_size, "color": MUT}
    return ax


def zone_patch(ax, hull, alpha=.16):
    import cartopy.crs as ccrs
    from matplotlib.patches import Polygon
    ax.add_patch(Polygon(hull, closed=True, facecolor=HOT, alpha=alpha,
                         edgecolor=HOT, linewidth=1.3, zorder=3,
                         transform=ccrs.PlateCarree()))


def fig_explained(day, hull, out, dpi):
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    met, lat, lon, saa = load(day)
    hrs = (met - met[0]) / 3600
    slices = orbit_slices(lat)

    fig = plt.figure(figsize=(11, 11.9), dpi=dpi, facecolor="white")
    gs = fig.add_gridspec(3, 1, height_ratios=[1, 1, .44], hspace=.46,
                          left=.07, right=.97, top=.94, bottom=.065)

    ax = basemap(fig.add_subplot(gs[0], projection=ccrs.PlateCarree()),
                 [-150, 80, -50, 40], 8)
    zone_patch(ax, hull, .30)
    sl = slices[0]
    for part in split_seam(lon[sl]):
        if len(part) > 1:
            ax.plot(lon[sl][part], lat[sl][part], color="#2f6fd0", lw=2,
                    transform=ccrs.PlateCarree(), zorder=4)
    ax.plot([], [], color="#2f6fd0", lw=2, label="one orbit (95 min)")
    ax.add_patch(Polygon([[0, 0]], closed=True, facecolor=HOT, alpha=.30,
                         edgecolor=HOT, label="the SAA"))
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    ax.set_title("A.   The SAA is a fixed region of the Earth. Fermi circles every "
                 "95 minutes,\n      never straying more than 26° from the "
                 "equator — straight through it.",
                 fontsize=11.5, color=INK, loc="left", pad=10)

    ax = basemap(fig.add_subplot(gs[1], projection=ccrs.PlateCarree()),
                 [-115, 40, -40, 14], 8)
    zone_patch(ax, hull, .20)
    cmap = plt.get_cmap("viridis")
    mins = []
    for k in range(5):
        sl = slices[k]
        lo_, la_, sa_ = lon[sl], lat[sl], saa[sl]
        colour = cmap(k / 4 * .85)
        for part in split_seam(lo_):
            if len(part) > 1:
                ax.plot(lo_[part], la_[part], color=colour, lw=1.5, alpha=.9,
                        transform=ccrs.PlateCarree(), zorder=4)
        if sa_.any():
            ax.plot(lo_[sa_], la_[sa_], color=HOT, lw=5, alpha=.95,
                    transform=ccrs.PlateCarree(), zorder=6,
                    solid_capstyle="round")
            j = int(np.argmax(lo_[sa_]))
            mins.append(sa_.sum() / 60)
            ax.annotate(f"orbit {k+1}", (lo_[sa_][j] + 2.5, la_[sa_][j] + 1.5),
                        transform=ccrs.PlateCarree(), ha="left", va="bottom",
                        fontsize=9.5, fontweight="bold", color=colour, zorder=8,
                        bbox=dict(boxstyle="round,pad=.18", fc="white",
                                  ec="none", alpha=.85))
    ax.plot([], [], color=HOT, lw=5, label="inside the SAA")
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    seq = ", ".join(f"{m:.0f}" for m in mins)
    ax.set_title("B.   Earth turns 23.8° under each orbit, so successive tracks "
                 f"step WEST\n      across the zone — and the crossings shorten: "
                 f"{seq} min.", fontsize=11.5, color=INK, loc="left", pad=10)

    ax = fig.add_subplot(gs[2])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(RULE)
    ax.axhspan(0, 1, color="#eef1f6", zorder=1)
    edge = np.diff(saa.astype(np.int8))
    entries = np.flatnonzero(edge == 1) + 1
    exits = np.flatnonzero(edge == -1) + 1
    n = 0
    for a in entries:
        later = exits[exits > a]
        if not later.size:
            continue
        b = later[0]
        n += 1
        ax.axvspan(hrs[a], hrs[b], color=HOT, alpha=.85, zorder=2)
        ax.text((hrs[a] + hrs[b]) / 2, .5, f"{(hrs[b]-hrs[a])*60:.0f}",
                ha="center", va="center", fontsize=7.5, color="white",
                fontweight="bold", zorder=3)
    gap_a = hrs[exits[5]] if exits.size > 5 else 9
    gap_b = hrs[entries[6]] if entries.size > 6 else 19
    ax.annotate("", xy=(gap_a, 1.14), xytext=(gap_b, 1.14),
                arrowprops=dict(arrowstyle="<->", color=MUT, lw=1))
    ax.text((gap_a + gap_b) / 2, 1.24,
            f"{gap_b-gap_a:.1f} h clear — the track has walked off the "
            f"western edge", ha="center", va="bottom", fontsize=8.5, color=MUT)
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xticks(range(0, 25, 3))
    ax.tick_params(colors=MUT, labelsize=8.5)
    ax.set_xlabel(f"hours UTC on 20{day[:2]}-{day[2:4]}-{day[4:]}",
                  fontsize=9.5, color=MUT)
    before = sum(1 for a in entries if hrs[a] < 10)
    ax.set_title(f"C.   {before} crossings in the first 9 hours, a "
                 f"{gap_b-gap_a:.0f}-hour clear stretch, then {n-before} more "
                 f"— {n} in all.\n      Numbers are minutes spent inside.",
                 fontsize=11.5, color=INK, loc="left", pad=26)
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return n


def fig_all_zoom(day, hull, out, dpi):
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    met, lat, lon, saa = load(day)
    slices = orbit_slices(lat)
    cmap = plt.get_cmap("turbo")
    extent = [-118, 42, -40, 14]

    fig = plt.figure(figsize=(13, 8), dpi=dpi, facecolor="white")
    ax = basemap(fig.add_axes([.055, .075, .875, .80],
                              projection=ccrs.PlateCarree()), extent)
    zone_patch(ax, hull, .13)
    spots = []
    for k, sl in enumerate(slices):
        lo_, la_, sa_ = lon[sl], lat[sl], saa[sl]
        colour = cmap(.04 + k / (len(slices) - 1) * .92)
        hit = sa_.any()
        for part in split_seam(lo_):
            if len(part) > 1:
                ax.plot(lo_[part], la_[part], color=colour,
                        lw=1.5 if hit else .9, alpha=.9 if hit else .30,
                        ls="-" if hit else (0, (3, 3)),
                        transform=ccrs.PlateCarree(), zorder=5 if hit else 4)
        if hit:
            # Each crossing in its own orbit's colour. A single red for all
            # ten is an unreadable crosshatch.
            for part in split_seam(lo_[sa_]):
                if len(part) > 1:
                    ax.plot(lo_[sa_][part], la_[sa_][part], color=colour, lw=6,
                            alpha=.95, transform=ccrs.PlateCarree(), zorder=6,
                            solid_capstyle="round",
                            path_effects=[pe.Stroke(linewidth=8.5,
                                                    foreground="white"),
                                          pe.Normal()])
            j = len(lo_[sa_]) // 2
            spots.append((float(lo_[sa_][j]), float(la_[sa_][j]), k + 1,
                          colour, sa_.sum() / 60))
    for x, y, k, colour, _ in spots:
        x = min(max(x, extent[0] + 6), extent[1] - 10)
        y = min(max(y, extent[2] + 4), extent[3] - 3)
        ax.annotate(f"{k}", (x, y), transform=ccrs.PlateCarree(), ha="center",
                    va="center", fontsize=10, fontweight="bold", color="white",
                    zorder=10,
                    path_effects=[pe.withStroke(linewidth=2.6, foreground=colour)])
    sm = ScalarMappable(Normalize(1, len(slices)), cmap)
    cb = fig.colorbar(sm, ax=ax, fraction=.026, pad=.012)
    cb.set_label("orbit number through the day", fontsize=9, color=MUT)
    cb.ax.tick_params(labelsize=8.5, colors=MUT)
    ax.plot([], [], color=MUT, lw=6,
            label="inside the SAA (thick, in that orbit's colour)")
    ax.plot([], [], color=MUT, lw=.9, ls=(0, (3, 3)), alpha=.5,
            label="orbit that misses the zone")
    ax.legend(loc="lower left", fontsize=9.5, frameon=False)
    n = int((np.diff(saa.astype(np.int8)) == -1).sum())
    order = "  ".join(f"{k}:{m:.0f}'" for _, _, k, _, m in
                      sorted(spots, key=lambda z: z[2]))
    ax.set_title(f"All {len(slices)} orbits of 20{day[:2]}-{day[2:4]}-{day[4:]} "
                 f"over the South Atlantic Anomaly\n{n} crossings. Numbers are "
                 f"the orbit; minutes inside: {order}",
                 fontsize=12, color=INK, loc="left", pad=12)
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return len(slices), n


def fig_all_small(day, hull, out, dpi):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon

    met, lat, lon, saa = load(day)
    slices = orbit_slices(lat)
    cmap = plt.get_cmap("turbo")
    rows = (len(slices) + 4) // 5
    fig, axes = plt.subplots(rows, 5, figsize=(12.5, 2.05 * rows), dpi=dpi,
                             facecolor="white")
    axes = np.atleast_2d(axes)
    for k in range(rows * 5):
        ax = axes[k // 5, k % 5]
        if k >= len(slices):
            ax.axis("off")
            continue
        ax.set_xlim(-180, 180)
        ax.set_ylim(-40, 40)
        ax.set_xticks([])
        ax.set_yticks([])
        for side in ax.spines.values():
            side.set_color("#dfe5ee")
        ax.add_patch(Polygon(hull, closed=True, facecolor=HOT, alpha=.18,
                             edgecolor=HOT, lw=.8))
        sl = slices[k]
        lo_, la_, sa_ = lon[sl], lat[sl], saa[sl]
        colour = cmap(.04 + k / (len(slices) - 1) * .92)
        for part in split_seam(lo_):
            if len(part) > 1:
                ax.plot(lo_[part], la_[part], color=colour, lw=1.3)
        mins = sa_.sum() / 60
        if sa_.any():
            for part in split_seam(lo_[sa_]):
                if len(part) > 1:
                    ax.plot(lo_[sa_][part], la_[sa_][part], color=colour, lw=3.4,
                            solid_capstyle="round")
        ax.set_title(f"orbit {k+1}" + (f"  ·  {mins:.0f} min" if mins > 0
                                       else "  ·  clear"),
                     fontsize=8.5, color=(INK if mins > 0 else MUT), pad=3)
    fig.tight_layout()
    fig.savefig(out, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--day", default="220827")
    parser.add_argument("--dpi", type=int, default=125)
    args = parser.parse_args()

    import matplotlib
    matplotlib.use("Agg")

    os.makedirs(OUTDIR, exist_ok=True)
    hull = saa_hull()
    print(f"SAA hull: {len(hull)} vertices", flush=True)

    n = fig_explained(args.day, hull, f"{OUTDIR}/orbit_explained.png", args.dpi)
    print(f"  orbit_explained.png   ({n} crossings)", flush=True)
    orbits, n = fig_all_zoom(args.day, hull, f"{OUTDIR}/orbit_all_zoom.png", args.dpi)
    print(f"  orbit_all_zoom.png    ({orbits} orbits, {n} crossings)", flush=True)
    fig_all_small(args.day, hull, f"{OUTDIR}/orbit_all_small.png", args.dpi)
    print(f"  orbit_all_small.png", flush=True)

    for name in sorted(os.listdir(OUTDIR)):
        if name.startswith("orbit_"):
            size = os.path.getsize(os.path.join(OUTDIR, name))
            print(f"    {name:26s} {size/1024:7.0f} KB")


if __name__ == "__main__":
    main()
