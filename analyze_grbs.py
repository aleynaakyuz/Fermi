"""
Reproduce the GBM ground localization for a random sample of GRBs and compare
each result against the official GBM HEALPix map.

This generalises deneme/deneme.py, which was hand-tuned for one short burst
(bn200826187), to run unattended over an arbitrary sample. The two things that
were hard-coded there have to be derived per burst here:

  source / background windows
      Taken from the catalog's T90 interval rather than chosen by eye, which is
      essentially what RoboBA does automatically. The background windows sit
      either side of T90 with a gap that scales with burst duration, clipped to
      the time span trigdat actually covers.

  systematic model
      Chosen to match the official map's LMETHOD keyword: 'Interactive' maps to
      hitl_model (human-in-the-loop), anything else to robo_ba_model. Applying
      the wrong one would compare against a differently-shaped distribution.

Only bursts that have an official glg_healpix_all file are eligible, since
without one there is nothing to compare to. That restricts the pool to roughly
2018 onward.

Outputs, per burst:
    site/figures/lc_<trigger>.png       lightcurve with the windows marked
    site/figures/sky_<trigger>.png      my 90% contour vs the official one
    site/results.json                   numbers behind the comparison
"""
import argparse
import csv
import json
import os
import traceback

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord

from gdt.missions.fermi.gbm.trigdat import Trigdat
from gdt.missions.fermi.time import Time
from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.plot.sky import EquatorialPlot, get_lonlat
from gdt.core.plot.plot import SkyPoints
from gdt.missions.fermi.gbm.localization import (
    GbmHealPix, hitl_model, robo_ba_model,
)
from gdt.missions.fermi.gbm.localization.dol.legacy_dol import legacy_DoL
from astropy.io import fits

DATADIR = "data/gbm"
MANIFEST = "catalog/fermigbrst.csv"
SITEDIR = "site"
FIGDIR = os.path.join(SITEDIR, "figures")
LOC_ERANGE = (50.0, 300.0)      # standard GBM localization band, keV
SHORT_LONG_CUT = 2.0            # s, conventional T90 boundary
MIN_NET_FRACTION = 0.5          # the S/N window must keep this share of net counts

# Suffix keeping each window mode's figures separate, set from --window.
FIG_SUFFIX = ""
MAKE_PLOTS = True


def figname(kind, trigger):
    return os.path.join(FIGDIR, f"{kind}_{trigger}{FIG_SUFFIX}.png")


# ----------------------------------------------------------------- windowing
def choose_timescale(t90):
    """trigdat binning: fine enough to resolve the burst, coarse enough for S/N."""
    if t90 < SHORT_LONG_CUT:
        return 64
    if t90 < 10.0:
        return 256
    return 1024


def choose_windows(t90, t90_start, span):
    """Derive source and background windows from the catalog T90 interval.

    Returns (timescale_ms, src_window, [bg_window, ...]).
    """
    timescale = choose_timescale(t90)
    dt = timescale / 1000.0

    # Snap the T90 interval outward onto the trigdat bin grid.
    src_lo = np.floor(t90_start / dt) * dt
    src_hi = np.ceil((t90_start + t90) / dt) * dt
    if src_hi - src_lo < dt:
        src_hi = src_lo + dt

    lo_limit, hi_limit = span[0] + 2.0, span[1] - 2.0
    src_lo = max(src_lo, lo_limit)
    src_hi = min(src_hi, hi_limit)

    # Keep background clear of the burst; widen the gap for long bursts so
    # tails do not leak into the background fit.
    gap = max(5.0, 0.25 * t90)
    pre_hi = src_lo - gap
    pre_lo = max(pre_hi - 60.0, lo_limit)
    post_lo = src_hi + gap
    post_hi = min(post_lo + 80.0, hi_limit)

    windows = []
    if pre_hi - pre_lo >= 15.0:
        windows.append((pre_lo, pre_hi))
    if post_hi - post_lo >= 15.0:
        windows.append((post_lo, post_hi))
    if not windows:
        raise ValueError("no room for a background window inside trigdat span")

    return timescale, (float(src_lo), float(src_hi)), windows


def peak_window(trigdat, timescale, t90_span, bg_times):
    """Pick the contiguous sub-interval with the best signal-to-noise.

    This is what a Burst Advocate does by eye: localize on the bright part of
    the burst rather than the whole T90. On a long, soft-tailed burst the tail
    contributes background without much signal, which both widens the
    statistical region and lets the centroid drift.

    Maximises sum(net) / sqrt(sum(background)) over every contiguous bin range
    inside the burst, on the summed triggered detectors in the localization
    band. Returns None if no interval has positive signal.
    """
    summed = trigdat.sum_detectors(trigdat.triggered_detectors,
                                   timescale=timescale).slice_energy(*(LOC_ERANGE,))
    fitter = BackgroundFitter.from_phaii(summed, Polynomial, time_ranges=bg_times)
    fitter.fit(order=1)
    background = fitter.interpolate_bins(summed.data.tstart, summed.data.tstop)

    counts = np.asarray(summed.data.counts).sum(axis=1).astype(float)
    bkg = np.asarray(background.counts).sum(axis=1).astype(float)
    tstart = np.asarray(summed.data.tstart, dtype=float)
    tstop = np.asarray(summed.data.tstop, dtype=float)

    # Search a little beyond T90 so a burst whose catalog interval is slightly
    # mis-centred can still be captured.
    margin = 0.15 * max(t90_span[1] - t90_span[0], 1.0)
    inside = np.where((tstop > t90_span[0] - margin) &
                      (tstart < t90_span[1] + margin))[0]
    if inside.size == 0:
        return None

    net = counts - bkg
    cum_net = np.concatenate([[0.0], np.cumsum(net)])
    cum_bkg = np.concatenate([[0.0], np.cumsum(np.clip(bkg, 0.0, None))])

    lo_idx, hi_idx = int(inside[0]), int(inside[-1]) + 1

    # Floor on how much of the burst may be thrown away, as insurance against
    # the window collapsing onto a sliver.
    #
    # Note this does NOT explain the short-burst regressions. bn251220985
    # (T90 0.38 s) trims to 0.19 s and still retains 79% of the net counts, yet
    # its localization degrades from 3.08 to 12.52 deg -- so the damage is not
    # lost signal. For short bursts the T90 interval is already the bright
    # interval (trimming buys little) while the shorter exposure leaves the
    # chi2 surface noisier, which is the more likely cause. Whether short
    # bursts should be trimmed at all is a question for the full sample, not
    # something to patch blind on three bursts.
    #
    # The constraint is kept because it costs nothing and guards the long-burst
    # case, and peak_net_fraction is recorded so the effect is measurable.
    total_net = max(cum_net[hi_idx] - cum_net[lo_idx], 0.0)
    required = MIN_NET_FRACTION * total_net if total_net > 0 else 0.0

    def search(min_signal):
        best_snr, best = -np.inf, None
        for i in range(lo_idx, hi_idx):
            for j in range(i + 1, hi_idx + 1):
                signal = cum_net[j] - cum_net[i]
                noise = cum_bkg[j] - cum_bkg[i]
                if signal <= 0.0 or noise <= 0.0 or signal < min_signal:
                    continue
                snr = signal / np.sqrt(noise)
                if snr > best_snr:
                    best_snr, best = snr, (i, j)
        return best_snr, best

    best_snr, best = search(required)
    if best is None:                       # constraint unsatisfiable, relax it
        best_snr, best = search(0.0)
    if best is None:
        return None

    i, j = best
    kept = (cum_net[j] - cum_net[i]) / total_net if total_net > 0 else float("nan")
    return float(tstart[i]), float(tstop[j - 1]), float(best_snr), float(kept)


# ------------------------------------------------------------------ analysis
def localize(trigger, row, window_mode="t90", artifacts=None):
    """Run the DoL fit for one burst and compare with the official map."""
    bdir = os.path.join(DATADIR, trigger)
    trigdat_file = next(
        f for f in sorted(os.listdir(bdir)) if f.startswith("glg_trigdat_")
    )
    healpix_file = next(
        f for f in sorted(os.listdir(bdir)) if f.startswith("glg_healpix_")
    )

    trigdat = Trigdat.open(os.path.join(bdir, trigdat_file))

    t90 = float(row["t90"])
    t90_start = float(row["t90_start"])

    probe = trigdat.to_phaii(trigdat._detectors[0], timescale=64)
    span = (float(probe.data.tstart[0]), float(probe.data.tstop[-1]))
    timescale, src_time, bg_times = choose_windows(t90, t90_start, span)

    # The background windows stay anchored to T90 either way; only the source
    # window changes, so the two modes differ in exactly one choice.
    t90_window, peak_snr, peak_kept = src_time, None, None
    if window_mode == "peak":
        found = peak_window(trigdat, timescale, src_time, bg_times)
        if found is not None:
            src_time = (found[0], found[1])
            peak_snr, peak_kept = found[2], found[3]

    # --- per-detector source and background counts (deneme.py step 2) -------
    # A detector that is disabled for a burst reads identically zero. The
    # polynomial fitter cannot fit that (it leaves its model unbound and raises
    # UnboundLocalError), so those detectors get zero background directly.
    # They still have to occupy their slot: legacy_DoL takes no detector mask
    # and expects all 14 in order.
    src_counts, src_exposure, bg_counts, bg_exposure, bg_trigdet = [], [], [], [], []
    dead_detectors = []
    for det in trigdat._detectors:
        phaii = trigdat.to_phaii(det, timescale=timescale)

        src_bin = phaii.data.integrate_time(*src_time)
        src_counts.append(src_bin.counts.astype(np.int32))
        src_exposure.append(src_bin.exposure)

        if phaii.data.counts.sum() == 0:
            dead_detectors.append(det)
            bg_counts.append(np.zeros_like(src_counts[-1]))
            bg_exposure.append(src_bin.exposure)
            continue

        fitter = BackgroundFitter.from_phaii(phaii, Polynomial, time_ranges=bg_times)
        fitter.fit(order=1)
        bg_rates = fitter.interpolate_bins(phaii.data.tstart, phaii.data.tstop)
        bg_bin = bg_rates.integrate_time(*src_time)
        bg_counts.append(bg_bin.counts.astype(np.int32))
        bg_exposure.append(bg_bin.exposure)

        if det in trigdat.triggered_detectors:
            bg_trigdet.append(bg_rates.slice_energy(*LOC_ERANGE))

    # --- energy bins and spacecraft state (step 3) --------------------------
    phaii = trigdat.to_phaii(trigdat._detectors[0], timescale=timescale)
    energies = np.concatenate([phaii.data.emin, [phaii.data.emax[-1]]])
    crange = [int(np.digitize(e, energies, right=True) - 1) for e in LOC_ERANGE]

    tcenter = trigdat.trigtime + 0.5 * sum(src_time)
    frame = trigdat.poshist.at(Time(tcenter, format="fermi"))
    scpos = frame.obsgeoloc.xyz.to_value("km")
    quaternion = frame.quaternion.scalar_last

    flight_ra = trigdat.headers["PRIMARY"]["RA_OBJ"]
    flight_dec = trigdat.headers["PRIMARY"]["DEC_OBJ"]

    # --- the DoL fit itself (step 4) ----------------------------------------
    # legacy_DoL wants RATES, not counts: deadtime_correct() internally does
    # `mrates = int32(sduration * mrates)` to recover counts. deneme.py passed
    # counts, which survives only because a common factor on both source and
    # background largely cancels in the fnorm normalisation -- for its 2 s
    # window the two conventions agree exactly. Over a 100 s window the same
    # mistake inflates chi2 past the internal 9999999 clamp, every grid point
    # ends up identical, and the search dies on an empty argmin.
    #
    # The detector mask is not ours to set: legacy_DoL(usedet=...) is ignored
    # by localize(), which recomputes it as `background rate > 100 counts/s`
    # over the 12 NaI. That already drops dead detectors.
    src_exp = float(np.mean(src_exposure))
    bg_exp = float(np.mean(bg_exposure))
    mrates = np.rint(np.array(src_counts) / src_exp).astype(np.int32)
    brates = np.rint(np.array(bg_counts) / bg_exp).astype(np.int32)

    dol = legacy_DoL()
    loc = dol.eval(
        crange, mrates, brates, src_exp, bg_exp,
        scpos, quaternion, energies, flight_ra, flight_dec, int(tcenter),
        scat_opt=1,
    )
    best_ra, best_dec = np.degrees([loc["best"]["ra"], loc["best"]["dec"]])
    az = float(np.degrees(loc["best"]["az"]))

    # --- systematic model matched to the official map's method (step 5) -----
    official_path = os.path.join(bdir, healpix_file)
    with fits.open(official_path) as hdul:
        head = hdul[0].header
        lmethod = str(head.get("LMETHOD", "unknown")).strip()
        off_err_rad = float(head.get("ERR_RAD", np.nan))
        loc_ver = str(head.get("LOC_VER", "")).strip()

    # convolve() takes the model function plus its arguments, not a called
    # model -- hitl_model(az) returns (sigmas, fracs), which is not callable.
    grb_type = "short" if t90 < SHORT_LONG_CUT else "long"
    if lmethod.lower().startswith("interactive"):
        model_func, model_args, model_name = hitl_model, (az,), "hitl_model"
    else:
        model_func, model_args, model_name = robo_ba_model, (grb_type,), "robo_ba_model"

    stat_map = dol.to_GbmHealPix(loc, frame)
    sys_map = stat_map.convolve(model_func, *model_args).remove_earth()

    official = GbmHealPix.open(official_path).remove_earth()

    mine = SkyCoord(*sys_map.centroid, unit="deg")
    off = SkyCoord(*official.centroid, unit="deg")
    cat = SkyCoord(float(row["ra"]), float(row["dec"]), unit="deg")

    result = {
        "trigger": trigger,
        "name": row["name"].strip(),
        "t90": t90,
        "t90_start": t90_start,
        "grb_type": grb_type,
        "fluence": _maybe_float(row.get("fluence")),
        "timescale_ms": timescale,
        "window_mode": window_mode,
        "t90_window": list(t90_window),
        "peak_snr": peak_snr,
        "peak_net_fraction": peak_kept,
        "src_duration": float(src_time[1] - src_time[0]),
        "src_time": list(src_time),
        "bg_times": [list(w) for w in bg_times],
        "crange": crange,
        "erange_kev": list(LOC_ERANGE),
        "triggered_detectors": list(trigdat.triggered_detectors),
        "dead_detectors": dead_detectors,
        "lmethod": lmethod,
        "loc_ver": loc_ver,
        "sys_model": model_name,
        "my_ra": float(mine.ra.deg), "my_dec": float(mine.dec.deg),
        "my_area90": float(sys_map.area(0.9)),
        "my_stat_area90": float(stat_map.area(0.9)),
        "dol_stat_err": float(loc["best"]["err"]),
        "official_ra": float(off.ra.deg), "official_dec": float(off.dec.deg),
        "official_area90": float(official.area(0.9)),
        "official_err_rad": off_err_rad,
        "catalog_ra": float(cat.ra.deg), "catalog_dec": float(cat.dec.deg),
        "catalog_err": _maybe_float(row.get("error_radius")),
        "flight_ra": float(flight_ra), "flight_dec": float(flight_dec),
        "sep_mine_official": float(mine.separation(off).deg),
        "sep_mine_catalog": float(mine.separation(cat).deg),
        "sep_official_catalog": float(off.separation(cat).deg),
        "official_in_my_map": float(sys_map.confidence(off.ra.deg, off.dec.deg)),
        "mine_in_official_map": float(official.confidence(mine.ra.deg, mine.dec.deg)),
    }

    if artifacts is not None:
        # Hand back the live objects so a caller can draw both window modes,
        # the official map and an external position on one set of axes.
        artifacts.update(sys_map=sys_map, stat_map=stat_map, official=official,
                         mine=mine, off=off, cat=cat, trigdat=trigdat,
                         timescale=timescale, src_time=src_time,
                         bg_times=bg_times, bg_trigdet=bg_trigdet)

    if MAKE_PLOTS:
        _plot_lightcurve(trigger, trigdat, timescale, src_time, bg_times,
                         bg_trigdet, result)
        _plot_skymap(trigger, sys_map, official, mine, off, result)
        _plot_zoom(trigger, sys_map, official, mine, off, cat, result)
    return result


def _maybe_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------- plots
def _plot_lightcurve(trigger, trigdat, timescale, src_time, bg_times, bg_trigdet, result):
    pad = max(10.0, 0.6 * (src_time[1] - src_time[0]))
    view = (max(bg_times[0][0] - 5, -140), min(bg_times[-1][1] + 5, 480))

    summed = trigdat.sum_detectors(trigdat.triggered_detectors, timescale=timescale)
    summed = summed.slice_energy(LOC_ERANGE).slice_time(view)

    fig = plt.figure(figsize=(8.2, 3.9), dpi=110)
    lcplot = Lightcurve(summed.to_lightcurve(), canvas=None)
    if bg_trigdet:
        try:
            lcplot.set_background(bg_trigdet[0].sum_time(bg_trigdet))
        except Exception:
            pass

    ax = plt.gca()
    for i, window in enumerate(bg_times):
        ax.axvspan(*window, color="#3a7fe8", alpha=0.18, lw=0,
                   label="background" if i == 0 else None)

    if result["window_mode"] == "peak":
        # Show what was discarded: the full T90 span outlined, the bright
        # sub-interval actually used for the fit filled in.
        t90_lo, t90_hi = result["t90_window"]
        ax.axvspan(t90_lo, t90_hi, facecolor="none", edgecolor="#b07d3a",
                   lw=1.0, ls="--", label="full T90 (not used)")
        ax.axvspan(*src_time, color="#e8833a", alpha=0.30, lw=0,
                   label="source (best S/N)")
    else:
        ax.axvspan(*src_time, color="#e8833a", alpha=0.22, lw=0,
                   label="source (T90)")

    mode_note = ("best-S/N window" if result["window_mode"] == "peak"
                 else "full T90 window")
    ax.set_title(
        f"{trigger}   T90 = {result['t90']:.2f} s ({result['grb_type']})   "
        f"{timescale} ms bins   {int(LOC_ERANGE[0])}-{int(LOC_ERANGE[1])} keV   "
        f"[{mode_note}]",
        fontsize=9,
    )
    ax.legend(fontsize=7, loc="upper right", framealpha=0.9)
    ax.tick_params(labelsize=8)
    ax.xaxis.label.set_size(9)
    ax.yaxis.label.set_size(9)
    plt.tight_layout()
    plt.savefig(figname("lc", trigger), bbox_inches="tight")
    plt.close("all")


OFFICIAL_COLOR = "#2f6fd0"
MINE_COLOR = "#d9534f"


def _contour_elements(skyplot):
    try:
        return list(skyplot.loc_contours.to_list())
    except Exception:
        return []


def _plot_skymap(trigger, sys_map, official, mine, off, result):
    """Both 90% contours on one all-sky frame, so the overlap is visible."""
    skyplot = EquatorialPlot()      # sets its own figsize internally

    # add_localization takes no colour argument, so each contour set is
    # recoloured right after it is added. Tracking artist identity is needed
    # because both calls register under the same '0.9_n' keys.
    skyplot.add_localization(official, clevels=[0.9], gradient=False,
                             sun=False, earth=True, detectors=[],
                             galactic_plane=True)
    first = _contour_elements(skyplot)
    for element in first:
        element.color = OFFICIAL_COLOR
    seen = {id(e) for e in first}

    skyplot.add_localization(sys_map, clevels=[0.9], gradient=False,
                             sun=False, earth=False, detectors=[],
                             galactic_plane=False)
    for element in _contour_elements(skyplot):
        if id(element) not in seen:
            element.color = MINE_COLOR

    _mark_centroid(skyplot, off, OFFICIAL_COLOR, "P")
    _mark_centroid(skyplot, mine, MINE_COLOR, "X")

    ax = skyplot.ax
    ax.set_title(
        f"{trigger}      red = this DoL reproduction      "
        f"blue = official GBM map      "
        f"centroids {result['sep_mine_official']:.2f}° apart",
        fontsize=9, pad=24,
    )
    plt.savefig(figname("sky", trigger), bbox_inches="tight")
    plt.close("all")


def _wrap_near(values, center):
    """Unwrap RA so a contour straddling 0/360 plots as one piece."""
    return center + (np.asarray(values) - center + 180.0) % 360.0 - 180.0


def _plot_zoom(trigger, sys_map, official, mine, off, cat, result):
    """Close-up of the two probability regions, where agreement is readable.

    The all-sky frame shows Earth occultation and overall geometry but squeezes
    the contours to near-invisibility; this panel is what actually shows how
    well the reproduction matches.
    """
    fig, ax = plt.subplots(figsize=(5.6, 5.0), dpi=110)
    center_ra, center_dec = off.ra.deg, off.dec.deg

    spans = []
    for hpx, color, label in ((official, OFFICIAL_COLOR, "official GBM"),
                              (sys_map, MINE_COLOR, "this reproduction")):
        for clevel, lw, ls in ((0.9, 1.8, "-"), (0.5, 1.2, "--")):
            try:
                paths = hpx.confidence_region_path(clevel)
            except Exception:
                continue
            for k, path in enumerate(paths):
                ra = _wrap_near(path[:, 0], center_ra)
                dec = path[:, 1]
                spans.append((ra, dec))
                ax.plot(
                    ra, dec, color=color, lw=lw, ls=ls,
                    label=(f"{label} 90%" if (clevel == 0.9 and k == 0)
                           else (f"{label} 50%" if (clevel == 0.5 and k == 0) else None)),
                )

    ax.plot(_wrap_near([off.ra.deg], center_ra), [off.dec.deg], marker="P",
            color=OFFICIAL_COLOR, ms=11, mec="white", mew=1.0, ls="none",
            label="official centroid")
    ax.plot(_wrap_near([mine.ra.deg], center_ra), [mine.dec.deg], marker="X",
            color=MINE_COLOR, ms=11, mec="white", mew=1.0, ls="none",
            label="my centroid")
    ax.plot(_wrap_near([cat.ra.deg], center_ra), [cat.dec.deg], marker="*",
            color="#2b2b2b", ms=13, mec="white", mew=0.8, ls="none",
            label="catalog position")

    if spans:
        all_ra = np.concatenate([s[0] for s in spans])
        all_dec = np.concatenate([s[1] for s in spans])
        pad = 0.08 * max(np.ptp(all_ra), np.ptp(all_dec), 2.0)
        ax.set_xlim(all_ra.min() - pad, all_ra.max() + pad)
        ax.set_ylim(all_dec.min() - pad, all_dec.max() + pad)

    ax.invert_xaxis()                                   # RA increases leftward
    ax.set_aspect(1.0 / max(np.cos(np.radians(center_dec)), 0.05))
    ax.set_xlabel("RA (deg)", fontsize=9)
    ax.set_ylabel("Dec (deg)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(alpha=0.25, lw=0.6)
    ax.set_title(f"{trigger}  —  centroids {result['sep_mine_official']:.2f}° apart",
                 fontsize=9)
    ax.legend(fontsize=7, loc="best", framealpha=0.9)
    plt.tight_layout()
    plt.savefig(figname("zoom", trigger), bbox_inches="tight")
    plt.close("all")


def _mark_centroid(skyplot, coord, color, marker):
    """Mark a centroid, projected exactly the way the contours are.

    add_localization wraps contour vertices as GCRS before transforming into
    the plot frame, so the centroids go through the same path to keep marker
    and contour in register.
    """
    gcrs = SkyCoord(coord.ra.deg, coord.dec.deg, frame="gcrs", unit="deg")
    lon, lat = get_lonlat(gcrs.transform_to(skyplot._astropy_frame))
    SkyPoints(
        lon.value, lat.value, skyplot.ax,
        flipped=skyplot._flipped, frame=skyplot._frame,
        color=color, marker=marker, s=55, edgecolor="white", linewidth=0.6,
        zorder=12,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=10, help="how many GRBs to sample")
    parser.add_argument("--seed", type=int, default=20261022, help="RNG seed")
    parser.add_argument("--triggers", type=str, default="",
                        help="explicit comma-separated triggers, bypasses sampling")
    parser.add_argument("--window", choices=("t90", "peak"), default="t90",
                        help="source window: whole T90, or the best-S/N sub-interval")
    parser.add_argument("--erange", type=str, default="50,300",
                        help="localization band in keV, e.g. 100,300 (trigdat has "
                             "8 coarse channels, so only a few distinct bands exist)")
    parser.add_argument("--noplots", action="store_true",
                        help="skip figure generation, for parameter sweeps")
    parser.add_argument("--out", type=str, default="",
                        help="results filename, defaults to results[_<mode>].json")
    args = parser.parse_args()

    global FIG_SUFFIX, LOC_ERANGE, MAKE_PLOTS
    LOC_ERANGE = tuple(float(v) for v in args.erange.split(","))
    MAKE_PLOTS = not args.noplots
    FIG_SUFFIX = "" if args.window == "t90" else f"_{args.window}"
    os.makedirs(FIGDIR, exist_ok=True)

    rows = {r["trigger_name"].strip(): r
            for r in csv.DictReader(open(MANIFEST, newline=""))}

    # Eligible = has an official healpix map and a usable T90.
    eligible = []
    for trigger, row in rows.items():
        bdir = os.path.join(DATADIR, trigger)
        if not os.path.isdir(bdir):
            continue
        if not any(f.startswith("glg_healpix_") for f in os.listdir(bdir)):
            continue
        if not any(f.startswith("glg_trigdat_") for f in os.listdir(bdir)):
            continue
        t90, t90_start = _maybe_float(row["t90"]), _maybe_float(row["t90_start"])
        if t90 is None or t90_start is None or not np.isfinite(t90) or t90 <= 0:
            continue
        eligible.append(trigger)
    eligible.sort()
    print(f"{len(eligible)} bursts eligible (official healpix + usable T90)")

    if args.triggers:
        sample = [t.strip() for t in args.triggers.split(",") if t.strip()]
    else:
        rng = np.random.default_rng(args.seed)
        sample = sorted(rng.choice(eligible, size=args.n, replace=False).tolist())
    print(f"sample ({len(sample)}): {', '.join(sample)}\n")

    results, failures = [], []
    for i, trigger in enumerate(sample, 1):
        print(f"[{i}/{len(sample)}] {trigger} ... ", end="", flush=True)
        try:
            result = localize(trigger, rows[trigger], window_mode=args.window)
            results.append(result)
            print(
                f"sep={result['sep_mine_official']:.2f} deg  "
                f"my_area={result['my_area90']:.0f}  "
                f"off_area={result['official_area90']:.0f} sq deg  "
                f"src={result['src_duration']:.1f}s"
            )
        except Exception as err:
            failures.append({"trigger": trigger, "error": f"{type(err).__name__}: {err}"})
            print(f"FAILED {type(err).__name__}: {err}")
            traceback.print_exc()

    payload = {
        "seed": args.seed,
        "window_mode": args.window,
        "n_eligible": len(eligible),
        "sample": sample,
        "results": results,
        "failures": failures,
    }
    name = args.out or (
        "results.json" if args.window == "t90" else f"results_{args.window}.json")
    with open(os.path.join(SITEDIR, name), "w") as handle:
        json.dump(payload, handle, indent=2)

    print(f"\n{len(results)} succeeded, {len(failures)} failed")
    if results:
        seps = np.array([r["sep_mine_official"] for r in results])
        print(f"separation from official map: median {np.median(seps):.2f} deg, "
              f"min {seps.min():.2f}, max {seps.max():.2f}")
    print(f"wrote {SITEDIR}/{name}")


if __name__ == "__main__":
    main()
