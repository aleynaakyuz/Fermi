"""
Build the single-page burst viewer: pick a class, pick a GRB, see its plot.

Every one of the ~2,000 bursts is selectable. That is only possible because
the sky regions are stored as contour polygons and drawn as SVG in the
browser; the PNG gallery it replaces totalled 81 MB against a 16 MB page
budget, so it could only ever show a few hundred.

Output: site/viewer.html
"""
import base64
import json
import glob
import os
import statistics
import datetime

import numpy as np

SRC = "site/vec/bursts"
LCSRC = "site/vec/lc"
DETSRC = "site/vec/det"
CHSRC = "site/vec/chan"
OUT = "site/viewer.html"

# Must match extract_channels.py.
CHAN_EMIN, CHAN_EMAX = 8.0, 900.0
CHAN_BANDS = [8, 15, 25, 50, 100, 200, 350, 600, 900]

ALL_DETS = [f"n{d}" for d in "0123456789ab"] + ["b0", "b1"]

# Per-burst annotations, shown in the panel when that burst is selected.
# Only for bursts where something about the event itself explains the result.
BURST_NOTES = {
    "bn200415367": (
        "Exceptionally strong and spectrally hard, and does not fit DoL's "
        "spectral templates. 19.7% of its net counts fall above 300 keV, the "
        "92nd percentile for short bursts, so the standard 50-300 keV band "
        "samples it poorly and the fit lands far off. The counts themselves "
        "are sound. Catalogued as a GRB, but its adopted position carries a "
        "0.12 deg error radius, meaning it came from another instrument."
    ),
    "bn250223692": (
        "Background fit unreliable. The spacecraft was in an extreme particle "
        "environment: 30,094 counts/s across all channels against a catalog "
        "median of 2,489, with livetime swinging between 45% and 100% as it "
        "left the region. The deadtime-corrected background is therefore "
        "genuinely non-stationary, and no polynomial order recovers it -- "
        "orders 1 to 4 were tried, and order 4 diverged to -368,000 counts/s. "
        "Treat both positions here as unmeasured rather than wrong."
    ),
    "bn200306934": (
        "The full-T90 result here is a non-detection, not a bad position. "
        "T90 is 0.99 s but T50 only 0.10 s, so the burst is a tenth-of-a-second "
        "spike inside a one-second T90. Integrating all of it multiplies the "
        "background 8.5-fold while adding no signal, and significance collapses "
        "from 5.7 sigma in the best-S/N window to 0.76 sigma over the full T90. "
        "With no signal to fit, the region inflates to 23,223 sq deg, 56% of the "
        "sky, so its 43.9 deg separation is a centroid drawn from half the sky "
        "rather than a measurement."
    ),
}

# A polynomial background assumes roughly constant counting conditions. Above
# this livetime swing that assumption is broken and the background estimate is
# not trustworthy, whatever the fit residual looks like. 7 of 2049 bursts.
LIVETIME_SWING_LIMIT = 0.20

OFFICIAL_C = "#2f6fd0"
T90_C = "#b8690a"
PEAK_C = "#d9534f"
EXT_C = "#1f9e6e"


def pack_contours(paths):
    """Pack contour polygons as base64 uint16 pairs at 0.01 degree.

    The four sky regions were the page's largest cost at 8.55 MB, written as
    JSON text at 12.8 characters per point against the four bytes a point
    actually needs. The stored coordinates already sit exactly on a 0.01
    degree grid -- checked across all 666,145 points, maximum deviation
    4e-12 -- so this throws nothing away.

    RA*100 reaches 36000 and Dec*100 spans -9000 to 9000, so Dec carries a
    +90 degree offset and both fit unsigned. Little-endian, to match the
    DataView read in the browser.
    """
    out = []
    for path in paths:
        a = np.asarray(path, dtype=float)
        q = np.empty(a.size, dtype="<u2")
        q[0::2] = np.rint(a[0::2] * 100.0)
        q[1::2] = np.rint((a[1::2] + 90.0) * 100.0)
        out.append(base64.b64encode(q.tobytes()).decode("ascii"))
    return out


def pack_lightcurve(row):
    """Re-encode a lightcurve record to cost less in the page.

    The page carries every burst inline against a 16 MB cap, and the summed
    lightcurve was its single largest component at 2.19 KB per burst. Two
    changes recover most of that without touching what is drawn:

      e -> es   trigdat's bin edges are piecewise uniform, so they are stored
                as [start, width, count] runs. The widths carry float jitter
                (1.024 against 1.025), which is why there are around 26 runs
                per burst rather than the four its documented template
                implies, but that still beats 144 explicit edges.
      r, b -> base64 uint16 scaled to the burst's own maximum. Count rates
                reach a few thousand per second, so 65,536 levels is far
                finer than the one-count rounding the arrays already had.

    Stored alongside as rs/bs is that maximum, so the browser can undo it.
    """
    edges = [float(v) for v in row["e"]]
    runs = []
    i = 0
    while i < len(edges) - 1:
        width = round(edges[i + 1] - edges[i], 4)
        n = 1
        while (i + n + 1 < len(edges)
               and round(edges[i + n + 1] - edges[i + n], 4) == width):
            n += 1
        runs.append([round(edges[i], 4), width, n])
        i += n

    def pack(values):
        arr = np.asarray(values, dtype=float)
        top = float(arr.max()) if arr.size else 0.0
        if top <= 0:
            top = 1.0
        q = np.clip(np.rint(arr / top * 65535.0), 0, 65535).astype("<u2")
        return base64.b64encode(q.tobytes()).decode("ascii"), round(top, 2)

    r_b64, r_top = pack(row["r"])
    b_b64, b_top = pack(row["b"])
    return {
        "ts": row["ts"], "es": runs,
        "r": r_b64, "rs": r_top, "b": b_b64, "bs": b_top,
        "w90": row["w90"], "wsn": row["wsn"], "wbg": row["wbg"],
    }


def main():
    records = []
    for path in sorted(glob.glob(os.path.join(SRC, "*.json"))):
        with open(path) as handle:
            records.append(json.load(handle))
    if not records:
        raise SystemExit("no contour data; run extract_contours.py first")

    for rec in records:
        for key in ("O", "A", "P", "R"):
            if rec.get(key):
                rec[key] = pack_contours(rec[key])

    # Lightcurves live in their own files because they are produced by a
    # separate, much cheaper pass that does not need the DoL fit.
    lc = {}
    for path in glob.glob(os.path.join(LCSRC, "*.json")):
        with open(path) as handle:
            row = json.load(handle)
        lc[row["t"]] = row
    missing = 0
    for rec in records:
        row = lc.get(rec["t"])
        if row is None:
            missing += 1
            continue
        rec["L"] = pack_lightcurve(row)
        live_min, swing, allrate = row.get("lv", [1.0, 0.0, 0])
        if swing > LIVETIME_SWING_LIMIT:
            rec["warn"] = (
                f"Unstable livetime: it swings by {swing:.0%} over the run "
                f"(low of {live_min:.0%}), with {allrate:,} counts/s across all "
                f"channels. A polynomial background assumes steady counting "
                f"conditions, so the background estimate here -- and both "
                f"positions drawn from it -- should be treated as unreliable."
            )
    print(f"  lightcurves attached to {len(records) - missing}/{len(records)} bursts")

    # Per-detector curves from TTE. A separate pass again, and one that only
    # covers bursts whose 14 TTE files are on disk, so the panel is drawn only
    # where there is something to draw.
    det = {}
    for path in glob.glob(os.path.join(DETSRC, "*.json")):
        with open(path) as handle:
            row = json.load(handle)
        det[row["t"]] = row
    for rec in records:
        row = det.get(rec["t"])
        if row is None:
            continue
        rec["D"] = {k: row[k] for k in ("lo", "hi", "nb", "sn", "sb",
                                        "d", "net")}
        # 2,032 of 2,049 bursts have all fourteen detectors in the usual
        # order, so the list is only carried for the handful that do not --
        # HEASARC omits TTE for the odd detector on the odd trigger. Spelling
        # it out every time cost 145 KB of a page with little room to spare.
        if row["dets"] != ALL_DETS:
            rec["D"]["dets"] = row["dets"]
        if "trg" in row:
            rec["D"]["trg"] = row["trg"]
        if "dead" in row:
            rec["D"]["dead"] = row["dead"]
    print(f"  detector panels attached to {sum(1 for r in records if 'D' in r)}"
          f"/{len(records)} bursts")

    chan = {}
    for path in glob.glob(os.path.join(CHSRC, "*.json")):
        with open(path) as handle:
            row = json.load(handle)
        chan[row["t"]] = row
    for rec in records:
        row = chan.get(rec["t"])
        if row is None:
            continue
        rec["C"] = {k: row[k] for k in ("lo", "hi", "nt", "nc", "sg", "sr",
                                        "bd", "bs", "dets", "tot")}
    print(f"  channel panels attached to {sum(1 for r in records if 'C' in r)}"
          f"/{len(records)} bursts")

    # The display bands are the same for every burst, so they are stated once
    # in the caption rather than repeated 2,000 times in the data.
    all_dets = json.dumps(ALL_DETS)
    ch_bands = json.dumps(CHAN_BANDS)
    ch_emin, ch_emax = CHAN_EMIN, CHAN_EMAX
    any_det = next(iter(det.values()), None)
    nai_band = "8-900 keV"
    if any_det:
        band = any_det["eb"]["nai"]
        nai_band = f"{band[0]:g}-{band[1]:g} keV"

    for rec in records:
        note = BURST_NOTES.get(rec["t"])
        if note:
            rec["note"] = note
    print(f"  notes attached to {sum(1 for r in records if 'note' in r)} bursts")
    print(f"  livetime flag raised on {sum(1 for r in records if 'warn' in r)} bursts")

    longs = [r for r in records if r["k"] == "long"]
    shorts = [r for r in records if r["k"] == "short"]

    def summary(rs):
        a = [r["as"] for r in rs]
        b = [r["ps"] for r in rs]
        ext = [r for r in rs if "x" in r]
        return {
            "n": len(rs),
            "meda": round(statistics.median(a), 2),
            "medb": round(statistics.median(b), 2),
            "imp": sum(1 for x, y in zip(a, b) if y < x),
            "ext": len(ext),
            "trutho": (round(statistics.median([r["x"]["do"] for r in ext]), 2)
                       if ext else None),
            "truthp": (round(statistics.median([r["x"]["dp"] for r in ext]), 2)
                       if ext else None),
        }

    stats = {"long": summary(longs), "short": summary(shorts)}
    data = json.dumps(records, separators=(",", ":"))
    built = datetime.date.today().isoformat()
    sl, ss = stats["long"], stats["short"]

    html = f"""<title>GRB Localization Viewer</title>
<style>
  :root {{
    --paper:#f5f7fa; --surface:#ffffff; --surface-2:#eef1f6;
    --ink:#141b26; --ink-soft:#46525f; --muted:#6b7686; --rule:#d8dee8;
    --accent:{OFFICIAL_C}; --t90:{T90_C}; --peak:{PEAK_C}; --ext:{EXT_C};
    --robo:#7c5cd6;
    --good:#2e7d55; --crit:#bf3b36; --grid:#c9d2de; --warn-note:#b0761c;
    --shadow:0 1px 2px rgba(20,27,38,.06), 0 6px 18px rgba(20,27,38,.05);
    --serif: ui-serif, Georgia, "Iowan Old Style", "Times New Roman", serif;
    --sans: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
    --mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
      --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
      --accent:#6fa3ea; --t90:#ffd166; --peak:#ef8079; --ext:#4cc79a;
      --robo:#b39aff;
      --good:#58b98a; --crit:#e8706a; --grid:#2f3b4d; --warn-note:#d8a04e;
      --shadow:0 1px 2px rgba(0,0,0,.4), 0 6px 18px rgba(0,0,0,.3);
    }}
  }}
  :root[data-theme="dark"] {{
    --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
    --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
    --accent:#6fa3ea; --t90:#ffd166; --peak:#ef8079; --ext:#4cc79a;
    --robo:#b39aff;
    --good:#58b98a; --crit:#e8706a; --grid:#2f3b4d; --warn-note:#d8a04e;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 6px 18px rgba(0,0,0,.3);
  }}
  *,*::before,*::after {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--paper); color:var(--ink);
    font-family:var(--sans); font-size:16px; line-height:1.6;
    -webkit-font-smoothing:antialiased; }}
  .wrap {{ max-width:1200px; margin:0 auto;
    padding:clamp(1.4rem,4vw,2.8rem) clamp(.8rem,3vw,1.6rem) 4rem;
    display:flex; flex-direction:column; gap:2rem; }}
  code {{ font-family:var(--mono); font-size:.9em; }}

  .eyebrow {{ font-family:var(--mono); font-size:.72rem; letter-spacing:.14em;
    text-transform:uppercase; color:var(--muted); margin:0; }}
  h1 {{ font-family:var(--serif); font-weight:600; margin:.45rem 0 0;
    font-size:clamp(1.6rem,4vw,2.3rem); line-height:1.14;
    letter-spacing:-.012em; text-wrap:balance; }}
  .standfirst {{ font-family:var(--serif); color:var(--ink-soft); max-width:62ch;
    font-size:clamp(1rem,2vw,1.1rem); margin:.8rem 0 0; }}
  h2 {{ font-family:var(--serif); font-size:1.26rem; font-weight:600; margin:0 0 .7rem; }}
  p {{ max-width:68ch; margin:0 0 .85rem; }}
  section p:last-child {{ margin-bottom:0; }}

  /* --------------------------------------------------------- picker */
  .picker {{ display:flex; flex-wrap:wrap; gap:1rem; align-items:flex-end;
    padding:1.1rem 1.2rem; background:var(--surface);
    border:1px solid var(--rule); border-radius:12px; box-shadow:var(--shadow); }}
  .field {{ display:flex; flex-direction:column; gap:.32rem; }}
  .field label {{ font-size:.72rem; letter-spacing:.08em; text-transform:uppercase;
    color:var(--muted); font-weight:650; }}
  select {{ font:inherit; font-size:.92rem; padding:.5rem .7rem;
    border:1px solid var(--rule); border-radius:8px;
    background:var(--paper); color:var(--ink); min-width:150px; }}
  select#grb {{ min-width:270px; font-family:var(--mono); font-size:.86rem; }}
  select:focus-visible, button:focus-visible {{
    outline:2px solid var(--accent); outline-offset:2px; }}
  .navbtns {{ display:flex; gap:.4rem; }}
  button.nav {{ font:inherit; font-size:.88rem; cursor:pointer;
    padding:.5rem .8rem; border-radius:8px; border:1px solid var(--rule);
    background:var(--paper); color:var(--ink); }}
  button.nav:hover {{ background:var(--surface-2); }}
  .pickcount {{ font-size:.8rem; color:var(--muted); margin-left:auto; }}

  /* ----------------------------------------------------------- plot */
  .stage {{ display:flex; flex-direction:column; gap:1.2rem; }}
  /* Plots get the full page width and sit side by side; the sky region is
     roughly square and the lightcurve wants length, so the split is uneven. */
  /* Two sky panels across the top, the lightcurve full width beneath them. */
  .plots {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr));
    gap:1.3rem; align-items:start; }}
  @media (max-width:860px) {{ .plots {{ grid-template-columns:1fr; }} }}
  .panelcap {{ font-size:.74rem; color:var(--muted); line-height:1.45;
    margin:.45rem 0 0; max-width:none; }}
  .plotbox {{ background:var(--surface); border:1px solid var(--rule);
    border-radius:12px; padding:1rem; box-shadow:var(--shadow); }}
  svg.sky {{ width:100%; height:auto; display:block; }}
  svg.sky text {{ font-family:var(--sans); fill:var(--muted); }}
  .axis line {{ stroke:var(--grid); stroke-width:.5; }}
  .frame {{ fill:none; stroke:var(--rule); stroke-width:1; }}

  .legend {{ display:flex; flex-wrap:wrap; gap:.75rem 1.1rem; margin-top:.8rem;
    font-size:.8rem; color:var(--ink-soft); }}
  .key {{ display:flex; align-items:center; gap:.4rem; }}
  .sw {{ width:20px; height:0; border-top-width:3px; border-top-style:solid; }}
  .sw.off {{ border-color:var(--accent); }}
  .sw.t90 {{ border-color:var(--t90); border-top-style:dashed; }}
  .sw.pk  {{ border-color:var(--peak); }}
  .sw.rb  {{ border-color:var(--robo); border-top-style:dashed; }}
  .star {{ color:var(--ext); font-size:1rem; line-height:1; }}
  .lcwrap {{ display:flex; flex-direction:column; margin-top:1rem;
    padding-top:.9rem; border-top:1px solid var(--rule); }}
  svg.lc {{ width:100%; height:auto; display:block; }}
  svg.lc text {{ font-family:var(--sans); fill:var(--muted); }}
  .band {{ width:18px; height:11px; border-radius:2px; }}
  /* Drawn with solid colour and a border rather than a faint tint: a 0.2
     alpha wash over the dark surface reads as mud. */
  .band.t90b {{ background:transparent; border:1.5px dashed var(--t90); }}
  .band.snb {{ background:var(--t90); }}
  .band.bgb {{ background:var(--accent); opacity:.45; }}

  /* ------------------------------------------------- per-detector panels */
  .detwrap {{ margin-top:1rem; padding-top:.9rem;
    border-top:1px solid var(--rule); }}
  /* Six columns, so the twelve NaIs fall as a tidy 6x2 and the two BGOs sit
     under them at exactly the same cell width. auto-fit was tried first and
     chose eight columns at this page width, which wraps the NaIs 8 + 4 and
     blows the BGO cells up to double size. The column count steps down by
     breakpoint instead. */
  .detgrid {{ display:grid; gap:.45rem; grid-template-columns:repeat(6,1fr); }}
  @media (max-width:980px) {{ .detgrid {{ grid-template-columns:repeat(4,1fr); }} }}
  @media (max-width:620px) {{ .detgrid {{ grid-template-columns:repeat(3,1fr); }} }}
  @media (max-width:420px) {{ .detgrid {{ grid-template-columns:repeat(2,1fr); }} }}
  /* Eight energy bands, so four columns rather than six: 6 + 2 leaves a
     ragged second row. */
  .detgrid.bands {{ grid-template-columns:repeat(4,1fr); }}
  @media (max-width:620px) {{ .detgrid.bands {{ grid-template-columns:repeat(2,1fr); }} }}
  svg.dsvg {{ width:100%; height:auto; display:block;
    background:var(--surface-2,transparent); border-radius:5px; }}
  svg.dsvg text {{ font-family:var(--mono); fill:var(--muted); }}
  svg.dsvg text.dlab {{ fill:var(--ink-soft); font-weight:650; }}
  svg.dsvg text.dlab.on {{ fill:var(--peak); }}
  svg.dsvg text.dnet {{ fill:var(--ink-soft); }}
  svg.dsvg text.dsc {{ fill:var(--muted); opacity:.75; }}
  .detsplit {{ font-size:.72rem; color:var(--muted); margin:.7rem 0 .4rem;
    padding-top:.5rem; border-top:1px dashed var(--rule); }}
  .detaxis {{ margin-top:.55rem; font-family:var(--mono); font-size:.7rem; }}
  svg.spec {{ width:100%; height:auto; display:block; margin-bottom:.2rem; }}
  svg.spec text {{ font-family:var(--sans); fill:var(--muted); }}

  /* ---------------------------------------------------------- facts */
  .facts {{ background:var(--surface); border:1px solid var(--rule);
    border-radius:12px; padding:1.1rem 1.2rem; box-shadow:var(--shadow); }}
  .facts h3 {{ margin:0 0 .2rem; font-size:1.05rem; font-family:var(--mono); }}
  .facts .sub {{ color:var(--muted); font-size:.82rem; margin:0 0 .9rem; }}
  dl.kv {{ margin:0; display:grid; grid-template-columns:auto 1fr; gap:.3rem .9rem; }}
  .factgrid {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(260px,1fr));
    gap:.4rem 2rem; align-items:start; }}
  .facthead {{ display:flex; flex-wrap:wrap; align-items:baseline;
    gap:.5rem 1.2rem; margin-bottom:.8rem; }}
  .facthead h3 {{ margin:0; }}
  dl.kv dt {{ color:var(--muted); font-size:.82rem; }}
  dl.kv dd {{ margin:0; font-size:.84rem; text-align:right;
    font-variant-numeric:tabular-nums; }}
  .sep {{ border-top:1px solid var(--rule); margin:.9rem 0; }}
  .big {{ font-family:var(--serif); font-size:1.5rem; line-height:1.1;
    font-variant-numeric:tabular-nums; }}
  .big.good {{ color:var(--good); }} .big.bad {{ color:var(--crit); }}
  .delta {{ font-size:.82rem; }}
  .delta.good {{ color:var(--good); }} .delta.bad {{ color:var(--crit); }}
  .facthead .big {{ margin-left:auto; }}
  .burstnote {{ margin:0 0 .9rem; padding:.7rem .9rem; font-size:.84rem;
    line-height:1.5; color:var(--ink-soft); background:var(--surface-2);
    border-left:3px solid var(--warn-note); border-radius:0 8px 8px 0; }}
  .burstnote strong {{ color:var(--ink); }}
  /* The flag is a stronger statement than a note: it says the measurement
     rests on an assumption that does not hold here. */
  .burstwarn {{ margin:0 0 .9rem; padding:.7rem .9rem; font-size:.84rem;
    line-height:1.5; color:var(--ink-soft); background:var(--surface-2);
    border-left:3px solid var(--crit); border-radius:0 8px 8px 0;
    max-width:none; }}
  .burstwarn strong {{ color:var(--crit); }}
  /* At full page width the external readings also run as columns, with the
     label above its value so the two stay visually paired. */
  .extbox.wide dl.kv {{ grid-template-columns:repeat(auto-fit,minmax(200px,1fr));
    gap:.2rem 2rem; }}
  .extbox.wide dl.kv dt {{ grid-row:1; }}
  .extbox.wide dl.kv dd {{ grid-row:2; text-align:left; font-size:.95rem; }}
  .extbox {{ background:var(--surface-2); border-left:3px solid var(--ext);
    border-radius:0 8px 8px 0; padding:.7rem .9rem; margin-top:.9rem;
    font-size:.84rem; }}
  .extbox strong {{ color:var(--ext); }}
  .noext {{ color:var(--muted); font-size:.82rem; margin-top:.9rem;
    font-style:italic; }}

  .tiles {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(160px,1fr));
    gap:1px; background:var(--rule); border:1px solid var(--rule);
    border-radius:10px; overflow:hidden; }}
  .tile {{ background:var(--surface); padding:.95rem 1.1rem; }}
  .tile .v {{ font-family:var(--serif); font-size:1.6rem; line-height:1.1;
    display:block; font-variant-numeric:tabular-nums; }}
  .tile .v .was {{ font-size:.92rem; color:var(--muted);
    text-decoration:line-through; margin-right:.3rem; }}
  .tile .k {{ font-size:.76rem; color:var(--muted); display:block; margin-top:.25rem; }}
  .callout {{ background:var(--surface-2); border-left:3px solid var(--accent);
    border-radius:0 8px 8px 0; padding:.95rem 1.15rem; }}
  .callout p:last-child {{ margin-bottom:0; }}
  footer {{ border-top:1px solid var(--rule); padding-top:1.2rem;
    color:var(--muted); font-size:.82rem; }}
  @media (prefers-reduced-motion: reduce) {{ * {{ animation:none!important; transition:none!important; }} }}
</style>

<div class="wrap">

  <header>
    <p class="eyebrow">Fermi GBM &middot; localization viewer</p>
    <h1>Every GBM burst, localized twice and checked against everyone else</h1>
    <p class="standfirst">
      Pick a burst class and a GRB. Each plot shows the official Fermi GBM
      region, this reproduction run on the full T90 interval, the same
      reproduction run on the brightest sub-interval, and &mdash; where another
      instrument caught the same burst &mdash; its arcsecond position, the one
      marker that does not depend on GBM at all.
    </p>
  </header>

  <div class="picker">
    <div class="field">
      <label for="cls">Burst class</label>
      <select id="cls">
        <option value="long">Long &mdash; T90 &ge; 2 s ({sl['n']})</option>
        <option value="short">Short &mdash; T90 &lt; 2 s ({ss['n']})</option>
      </select>
    </div>
    <div class="field">
      <label for="grb">GRB</label>
      <select id="grb"></select>
    </div>
    <div class="field">
      <label for="sort">Order by</label>
      <select id="sort">
        <option value="name">Name</option>
        <option value="best">Best agreement first</option>
        <option value="worst">Worst agreement first</option>
        <option value="gain">Biggest gain from trimming</option>
        <option value="ext">Has external position</option>
      </select>
    </div>
    <div class="navbtns">
      <button class="nav" id="prev" type="button">&larr; Prev</button>
      <button class="nav" id="next" type="button">Next &rarr;</button>
    </div>
    <span class="pickcount" id="pickcount"></span>
  </div>

  <div class="stage">
    <div class="plotbox">
     <div class="plots">
      <div>
      <svg class="sky" id="sky" viewBox="0 0 440 400" role="img"
           aria-label="Sky localization in right ascension and declination"></svg>
      <p class="panelcap">Right ascension and declination. Familiar
        coordinates, but stretched near the poles, where a region containing a
        pole cannot close.</p>
      </div>
      <div>
      <svg class="sky" id="sky2" viewBox="0 0 440 400" role="img"
           aria-label="Sky localization projected about the burst, showing true shape"></svg>
      <p class="panelcap">Projected about the burst. True shape and size at any
        declination; rings mark degrees from the centre.</p>
      </div>
     </div>

      <div class="legend">
        <span class="key"><span class="sw off"></span> official Fermi GBM 90%</span>
        <span class="key"><span class="sw t90"></span> mine, full T90</span>
        <span class="key"><span class="sw pk"></span> mine, best S/N</span>
        <span class="key"><span class="sw rb"></span> same fit, RoboBA error model</span>
        <span class="key"><span class="star">&#9733;</span> other mission</span>
      </div>
      <p class="panelcap" style="max-width:68ch">
        The purple region is <em>not</em> a separate localization. It is the
        best-S/N fit convolved with RoboBA's systematic kernel instead of the
        human-in-the-loop one. Both kernels are symmetric, so the centre barely
        moves and essentially only the size changes &mdash; a median 0.66 times
        the area, because RoboBA's kernel carries no 13&ndash;15&deg; tail. It
        is here to show how much of a region's size is the assumed error model
        rather than the data. RoboBA's own interval selection is not part of
        this, and is not public.
      </p>

      <div class="lcwrap">
        <svg class="lc" id="lc" viewBox="0 0 960 300" role="img"
             aria-label="Lightcurve for the selected burst with the source and background windows marked"></svg>
        <div class="legend">
          <span class="key"><span class="band t90b"></span> full T90</span>
          <span class="key"><span class="band snb"></span> best-S/N window</span>
          <span class="key"><span class="band bgb"></span> background</span>
          <span class="key"><span class="sw pk"></span> fitted background</span>
        </div>
      </div>

      <div class="detwrap">
        <p class="panelcap" style="max-width:78ch;margin-bottom:.7rem">
          Each detector on its own, from TTE. GBM has no optics &mdash; it
          works out where a burst was purely by comparing how much each
          detector saw, so this is the fit's actual input. All twelve NaI
          panels share one vertical scale, printed at the lower right of each,
          because the comparison <em>between</em> them is the whole signal;
          scaling each to its own peak would make every detector look alike.
          The two BGOs are scaled separately, since their band is far wider
          and their background correspondingly higher. Detectors that
          triggered are marked <span style="color:var(--peak)">&#9733;</span>,
          and the figure at the top right is net counts in the T90 window.
        </p>
        <div id="dets"></div>
      </div>

      <div class="detwrap">
        <p class="panelcap" style="max-width:78ch;margin-bottom:.7rem">
          The same burst resolved in energy. A GRB is not static spectrally:
          most decay hard to soft, pulses peak earlier at high energy than at
          low, and a hard precursor can be invisible in the summed curve. The
          map is shaded by significance above background,
          (counts &minus; background) / &radic;background, rather than by raw
          rate &mdash; the background falls by more than two orders of
          magnitude from 8 to 900 keV, so a raw map would just be a picture
          of the background with the burst lost in it.
        </p>
        <div id="chan"></div>
      </div>
    </div>

    <div class="facts" id="facts"></div>
  </div>

  <section>
    <h2>How the two classes compare</h2>
    <div class="tiles">
      <div class="tile"><span class="v"><span class="was">{sl['meda']:.2f}&deg;</span>{sl['medb']:.2f}&deg;</span>
        <span class="k">long: median separation from official</span></div>
      <div class="tile"><span class="v"><span class="was">{ss['meda']:.2f}&deg;</span>{ss['medb']:.2f}&deg;</span>
        <span class="k">short: median separation from official</span></div>
      <div class="tile"><span class="v">{sl['imp']}/{sl['n']}</span>
        <span class="k">long: improved by trimming</span></div>
      <div class="tile"><span class="v">{ss['imp']}/{ss['n']}</span>
        <span class="k">short: improved by trimming</span></div>
    </div>
    <p style="margin-top:1rem">
      Trimming the source window to the brightest part of the burst helps both
      classes, and helps long bursts more reliably &mdash; a short burst has
      almost no faint tail to discard, because its T90 interval already is the
      bright interval.
    </p>
  </section>

  <section>
    <h2>Measured against actual truth</h2>
    <div class="callout">
      <p>
        Separation from the official map measures agreement with GBM, not
        accuracy. For the {sl['ext'] + ss['ext']} bursts where another
        instrument supplied an arcsecond position, both can be scored against
        something independent: the official Fermi map sits a median
        {sl['trutho']:.2f}&deg; from truth for long bursts, this reproduction
        {sl['truthp']:.2f}&deg;. They are not meaningfully distinguishable in
        accuracy &mdash; the reproduction is as good at finding the burst as
        the pipeline it reproduces.
      </p>
    </div>
    <p style="margin-top:1rem">
      Select a burst with a star to see this directly: the panel reports how
      far the official map, the T90 window and the best-S/N window each landed
      from the external position. Agreeing with GBM and being right turn out to
      be different things, and occasionally the "worse" window is the one
      closer to the truth.
    </p>
  </section>

  <footer>
    Built {built} from {len(records)} bursts with an official GBM HEALPix map.
    Regions are 90% confidence, drawn from contour polygons rather than images
    so every burst is selectable. Localization by <code>legacy_DoL</code> from
    <code>astro-gdt-fermi</code>, 50&ndash;300&thinsp;keV, scattering on.
  </footer>
</div>

<script>
(function () {{
  var DATA = {data};
  var byClass = {{long: [], short: []}};
  DATA.forEach(function (r) {{ byClass[r.k].push(r); }});

  var clsSel = document.getElementById('cls');
  var grbSel = document.getElementById('grb');
  var sortSel = document.getElementById('sort');
  var svg = document.getElementById('sky');
  var facts = document.getElementById('facts');
  var pickcount = document.getElementById('pickcount');

  // Read at draw time, not load time, so the plot and the legend swatches stay
  // in agreement if the viewer's theme changes under us.
  function colors() {{
    var cs = getComputedStyle(document.documentElement);
    function v(name, fallback) {{
      return (cs.getPropertyValue(name) || '').trim() || fallback;
    }}
    return {{
      off: v('--accent', '{OFFICIAL_C}'),
      t90: v('--t90', '{T90_C}'),
      peak: v('--peak', '{PEAK_C}'),
      ext: v('--ext', '{EXT_C}'),
      robo: v('--robo', '#7c5cd6'),
      muted: v('--muted', '#6b7686')
    }};
  }}

  function current() {{ return byClass[clsSel.value]; }}

  function ordered() {{
    var list = current().slice();
    var how = sortSel.value;
    if (how === 'best') list.sort(function (a, b) {{ return a.ps - b.ps; }});
    else if (how === 'worst') list.sort(function (a, b) {{ return b.ps - a.ps; }});
    else if (how === 'gain') list.sort(function (a, b) {{ return (a.ps - a.as) - (b.ps - b.as); }});
    else if (how === 'ext') list.sort(function (a, b) {{
      var ax = a.x ? 0 : 1, bx = b.x ? 0 : 1;
      return ax - bx || a.t.localeCompare(b.t);
    }});
    else list.sort(function (a, b) {{ return a.t.localeCompare(b.t); }});
    return list;
  }}

  function fillGrbs(keep) {{
    var list = ordered();
    grbSel.innerHTML = list.map(function (r) {{
      var tag = (r.x ? ' *' : '') + (r.warn ? ' ~' : '') + (r.note ? ' !' : '');
      return '<option value="' + r.t + '">' + r.t + '  (' + r.ps.toFixed(2) +
             '\\u00B0' + tag + ')</option>';
    }}).join('');
    if (keep && list.some(function (r) {{ return r.t === keep; }})) grbSel.value = keep;
    pickcount.textContent = list.length +
      ' bursts  \\u00B7  * = external position  \\u00B7  ! = note  ' +
      '\\u00B7  ~ = unreliable background';
  }}

  // ---------------------------------------------------------------- plot
  function wrapNear(v, c) {{ return c + ((v - c + 180) % 360 + 360) % 360 - 180; }}

  function draw(r) {{
    unpackContours(r);
    var COL = colors();
    var W = 440, H = 400, M = {{l: 52, r: 14, t: 14, b: 40}};
    var cx = r.oc[0], cy = r.oc[1];
    var xs = [], ys = [];

    function collect(paths) {{
      paths.forEach(function (flat) {{
        for (var i = 0; i < flat.length; i += 2) {{
          xs.push(wrapNear(flat[i], cx));
          ys.push(flat[i + 1]);
        }}
      }});
    }}
    collect(r.O); collect(r.A); collect(r.P);
    [r.oc, r.ac, r.pc].forEach(function (c) {{ xs.push(wrapNear(c[0], cx)); ys.push(c[1]); }});
    if (r.x) {{ xs.push(wrapNear(r.x.c[0], cx)); ys.push(r.x.c[1]); }}
    if (!xs.length) {{ svg.innerHTML = ''; return; }}

    var x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
    var y0 = Math.min.apply(null, ys), y1 = Math.max.apply(null, ys);
    var padx = Math.max((x1 - x0) * 0.10, 0.6), pady = Math.max((y1 - y0) * 0.10, 0.6);
    x0 -= padx; x1 += padx; y0 -= pady; y1 += pady;

    var iw = W - M.l - M.r, ih = H - M.t - M.b;

    // Keep degrees square. One degree of RA covers cos(dec) times as much sky
    // as one degree of Dec, so without this the regions stretch to fill the
    // box and a circular error region reads as an ellipse.
    var cosd = Math.max(Math.cos(cy * Math.PI / 180), 0.05);
    var spanX = (x1 - x0) * cosd, spanY = y1 - y0;
    var want = iw / ih;
    if (spanX / spanY < want) {{
      var grow = (spanY * want - spanX) / cosd / 2;
      x0 -= grow; x1 += grow;
    }} else {{
      var growY = (spanX / want - spanY) / 2;
      y0 -= growY; y1 += growY;
    }}
    // RA increases to the left, the astronomical convention.
    function px(v) {{ return M.l + iw * (x1 - v) / (x1 - x0 || 1); }}
    function py(v) {{ return M.t + ih * (y1 - v) / (y1 - y0 || 1); }}

    var out = [];
    out.push('<rect class="frame" x="' + M.l + '" y="' + M.t +
             '" width="' + iw + '" height="' + ih + '" />');

    function ticks(lo, hi) {{
      var span = hi - lo, raw = span / 4;
      var mag = Math.pow(10, Math.floor(Math.log10(raw)));
      var step = [1, 2, 2.5, 5, 10].map(function (m) {{ return m * mag; }})
                  .filter(function (s) {{ return s >= raw; }})[0] || mag * 10;
      var t = [], v = Math.ceil(lo / step) * step;
      for (; v <= hi; v += step) t.push(Math.round(v * 1e6) / 1e6);
      return t;
    }}

    out.push('<g class="axis">');
    ticks(x0, x1).forEach(function (v) {{
      var X = px(v);
      out.push('<line x1="' + X + '" y1="' + M.t + '" x2="' + X + '" y2="' + (M.t + ih) + '"/>');
      out.push('<text x="' + X + '" y="' + (M.t + ih + 15) + '" font-size="9" text-anchor="middle">' +
               v.toFixed(Math.abs(v) < 10 ? 1 : 0) + '</text>');
    }});
    ticks(y0, y1).forEach(function (v) {{
      var Y = py(v);
      out.push('<line x1="' + M.l + '" y1="' + Y + '" x2="' + (M.l + iw) + '" y2="' + Y + '"/>');
      out.push('<text x="' + (M.l - 6) + '" y="' + (Y + 3) + '" font-size="9" text-anchor="end">' +
               v.toFixed(Math.abs(v) < 10 ? 1 : 0) + '</text>');
    }});
    out.push('</g>');
    out.push('<text x="' + (M.l + iw / 2) + '" y="' + (H - 8) +
             '" font-size="10" text-anchor="middle">RA (deg)</text>');
    out.push('<text transform="translate(13,' + (M.t + ih / 2) +
             ') rotate(-90)" font-size="10" text-anchor="middle">Dec (deg)</text>');

    function poly(paths, color, dash, width) {{
      paths.forEach(function (flat) {{
        var d = '';
        for (var i = 0; i < flat.length; i += 2) {{
          d += (i ? 'L' : 'M') + px(wrapNear(flat[i], cx)).toFixed(1) + ',' +
               py(flat[i + 1]).toFixed(1);
        }}
        if (!d) return;
        // Not closed with Z: a region crossing the seam comes back as several
        // pieces, and closing each one would draw a false edge along it. A
        // genuinely closed contour already repeats its first point.
        out.push('<path d="' + d + '" fill="none" stroke="' + color +
                 '" stroke-width="' + width + '"' +
                 (dash ? ' stroke-dasharray="5,3"' : '') +
                 ' stroke-linejoin="round" />');
      }});
    }}
    poly(r.O, COL.off, false, 2.1);
    poly(r.A, COL.t90, true, 2.0);
    poly(r.P, COL.peak, false, 1.6);
    if (r.R && r.R.length) poly(r.R, COL.robo, true, 1.4);

    function cross(c, color, kind) {{
      var X = px(wrapNear(c[0], cx)), Y = py(c[1]), s = 5;
      if (kind === 'plus') {{
        out.push('<path d="M' + (X - s) + ',' + Y + 'H' + (X + s) + 'M' + X + ',' + (Y - s) +
                 'V' + (Y + s) + '" stroke="' + color + '" stroke-width="2.2" />');
      }} else {{
        out.push('<path d="M' + (X - s) + ',' + (Y - s) + 'L' + (X + s) + ',' + (Y + s) +
                 'M' + (X + s) + ',' + (Y - s) + 'L' + (X - s) + ',' + (Y + s) +
                 '" stroke="' + color + '" stroke-width="2.2" />');
      }}
    }}
    cross(r.oc, COL.off, 'plus');
    cross(r.pc, COL.peak, 'x');

    if (r.x) {{
      var X = px(wrapNear(r.x.c[0], cx)), Y = py(r.x.c[1]);
      var pts = [], R = 8.5, rr = 3.6;
      for (var k = 0; k < 10; k++) {{
        var ang = -Math.PI / 2 + k * Math.PI / 5, rad = (k % 2 ? rr : R);
        pts.push((X + rad * Math.cos(ang)).toFixed(1) + ',' + (Y + rad * Math.sin(ang)).toFixed(1));
      }}
      out.push('<polygon points="' + pts.join(' ') + '" fill="' + COL.ext +
               '" stroke="#fff" stroke-width="1" />');
    }}
    svg.innerHTML = out.join('');
  }}

  // ----------------------------------------- burst-centred azimuthal plot
  // Azimuthal equidistant about the burst: distance from the centre is true
  // angular separation in any direction, so a region keeps its shape and size
  // wherever it sits. Unlike RA/Dec this has no singularity at the poles and
  // no seam at RA 0, which is what breaks the left-hand panel for the 3% of
  // bursts whose region reaches a pole.
  function drawAzim(r) {{
    unpackContours(r);
    var svg2 = document.getElementById('sky2');
    var COL = colors();
    var W = 440, H = 400, M = 26;
    var cx = W / 2, cy = H / 2;
    var ra0 = r.oc[0] * Math.PI / 180, dec0 = r.oc[1] * Math.PI / 180;
    var sin0 = Math.sin(dec0), cos0 = Math.cos(dec0);

    // (ra,dec) -> offset in degrees: north up, east left.
    function project(raDeg, decDeg) {{
      var ra = raDeg * Math.PI / 180, dec = decDeg * Math.PI / 180;
      var dra = ra - ra0;
      var cosc = sin0 * Math.sin(dec) + cos0 * Math.cos(dec) * Math.cos(dra);
      cosc = Math.max(-1, Math.min(1, cosc));
      var c = Math.acos(cosc) * 180 / Math.PI;
      var az = Math.atan2(Math.cos(dec) * Math.sin(dra),
                          cos0 * Math.sin(dec) - sin0 * Math.cos(dec) * Math.cos(dra));
      return [-c * Math.sin(az), c * Math.cos(az)];
    }}

    var maxr = 0, pts = [];
    function scan(paths) {{
      paths.forEach(function (flat) {{
        for (var i = 0; i < flat.length; i += 2) {{
          var q = project(flat[i], flat[i + 1]);
          pts.push(q);
          maxr = Math.max(maxr, Math.hypot(q[0], q[1]));
        }}
      }});
    }}
    scan(r.O); scan(r.A); scan(r.P);
    if (r.x) {{
      var qx = project(r.x.c[0], r.x.c[1]);
      maxr = Math.max(maxr, Math.hypot(qx[0], qx[1]));
    }}
    if (!(maxr > 0)) {{ svg2.innerHTML = ''; return; }}
    maxr *= 1.08;
    var k = (Math.min(W, H) / 2 - M) / maxr;        // pixels per degree
    function sx(q) {{ return cx + q[0] * k; }}
    function sy(q) {{ return cy - q[1] * k; }}

    var out = [];

    // rings at a round angular step, labelled in degrees from the centre
    var raw = maxr / 3;
    var mag = Math.pow(10, Math.floor(Math.log10(raw)));
    var step = [1, 2, 2.5, 5, 10].map(function (m) {{ return m * mag; }})
                .filter(function (v) {{ return v >= raw; }})[0] || mag * 10;
    out.push('<g class="axis">');
    for (var ring = step; ring <= maxr; ring += step) {{
      out.push('<circle cx="' + cx + '" cy="' + cy + '" r="' + (ring * k).toFixed(1) +
               '" fill="none" />');
      out.push('<text x="' + (cx + 3) + '" y="' + (cy - ring * k - 3).toFixed(1) +
               '" font-size="9">' + (ring % 1 ? ring.toFixed(1) : ring) +
               '\\u00B0</text>');
    }}
    for (var a = 0; a < 360; a += 45) {{
      var rad = a * Math.PI / 180;
      out.push('<line x1="' + cx + '" y1="' + cy +
               '" x2="' + (cx + Math.sin(rad) * maxr * k).toFixed(1) +
               '" y2="' + (cy - Math.cos(rad) * maxr * k).toFixed(1) + '"/>');
    }}
    out.push('</g>');
    out.push('<text x="' + cx + '" y="' + (M - 10) +
             '" font-size="10" text-anchor="middle">N</text>');
    out.push('<text x="' + (M - 12) + '" y="' + (cy + 4) +
             '" font-size="10" text-anchor="middle">E</text>');

    function poly(paths, color, dash, width) {{
      // A layer stored as several pieces is either split at the RA seam or
      // genuinely multi-island. The projection rejoins a seam split visually,
      // but closing each piece with Z would draw a false chord across it, so
      // only a single-piece layer is closed and filled.
      var solo = paths.length === 1;
      paths.forEach(function (flat) {{
        var d = '';
        for (var i = 0; i < flat.length; i += 2) {{
          var q = project(flat[i], flat[i + 1]);
          d += (i ? 'L' : 'M') + sx(q).toFixed(1) + ',' + sy(q).toFixed(1);
        }}
        if (!d) return;
        out.push('<path d="' + d + (solo ? 'Z' : '') + '" fill="' +
                 (solo ? color : 'none') + '"' +
                 (solo ? ' fill-opacity="0.07"' : '') +
                 ' stroke="' + color + '" stroke-width="' + width + '"' +
                 (dash ? ' stroke-dasharray="5,3"' : '') +
                 ' stroke-linejoin="round" />');
      }});
    }}
    poly(r.O, COL.off, false, 2.1);
    poly(r.A, COL.t90, true, 2.0);
    poly(r.P, COL.peak, false, 1.6);
    if (r.R && r.R.length) poly(r.R, COL.robo, true, 1.4);

    function mark(c, color, kind) {{
      var q = project(c[0], c[1]), X = sx(q), Y = sy(q), t = 5;
      if (kind === 'plus') {{
        out.push('<path d="M' + (X - t) + ',' + Y + 'H' + (X + t) + 'M' + X + ',' +
                 (Y - t) + 'V' + (Y + t) + '" stroke="' + color + '" stroke-width="2.2"/>');
      }} else {{
        out.push('<path d="M' + (X - t) + ',' + (Y - t) + 'L' + (X + t) + ',' + (Y + t) +
                 'M' + (X + t) + ',' + (Y - t) + 'L' + (X - t) + ',' + (Y + t) +
                 '" stroke="' + color + '" stroke-width="2.2"/>');
      }}
    }}
    mark(r.oc, COL.off, 'plus');
    mark(r.pc, COL.peak, 'x');

    if (r.x) {{
      var q = project(r.x.c[0], r.x.c[1]), X = sx(q), Y = sy(q);
      var star = [], R = 8.5, rr = 3.6;
      for (var n = 0; n < 10; n++) {{
        var ang = -Math.PI / 2 + n * Math.PI / 5, rad2 = (n % 2 ? rr : R);
        star.push((X + rad2 * Math.cos(ang)).toFixed(1) + ',' +
                  (Y + rad2 * Math.sin(ang)).toFixed(1));
      }}
      out.push('<polygon points="' + star.join(' ') + '" fill="' + COL.ext +
               '" stroke="#fff" stroke-width="1"/>');
    }}
    svg2.innerHTML = out.join('');
  }}

  // ------------------------------------------------- compact array decoding
  // The page carries every burst inline, so the bulky arrays travel encoded
  // and are expanded the first time a burst is actually looked at.
  function b64bytes(s) {{
    var bin = atob(s), n = bin.length, out = new Uint8Array(n);
    for (var i = 0; i < n; i++) out[i] = bin.charCodeAt(i);
    return out;
  }}

  // Contour polygons arrive as base64 uint16 pairs on a 0.01 degree grid:
  // RA directly, Dec offset by +90 so both fit unsigned. Expanded once per
  // burst and cached on the record, since both sky panels draw the same
  // polygons.
  function unpackContours(r) {{
    if (r.cu) return;
    ['O', 'A', 'P', 'R'].forEach(function (k) {{
      var src = r[k];
      if (!src || !src.length) return;
      r[k] = src.map(function (b64) {{
        var u = b64bytes(b64), n = u.length >> 1, out = new Float64Array(n);
        for (var i = 0; i < n; i++) {{
          var v = u[2 * i] | (u[2 * i + 1] << 8);
          out[i] = (i & 1) ? (v / 100 - 90) : (v / 100);
        }}
        return out;
      }});
    }});
    r.cu = 1;
  }}

  function unpackLC(L) {{
    if (L.e) return;                       // already expanded
    var e = [], i, j, seg;
    for (i = 0; i < L.es.length; i++) {{
      seg = L.es[i];
      for (j = 0; j < seg[2]; j++) e.push(seg[0] + seg[1] * j);
    }}
    seg = L.es[L.es.length - 1];
    e.push(seg[0] + seg[1] * seg[2]);
    L.e = e;
    function de(b64, top) {{
      var u = b64bytes(b64), n = u.length >> 1, out = new Array(n);
      for (var k = 0; k < n; k++) {{
        out[k] = (u[2 * k] | (u[2 * k + 1] << 8)) / 65535 * top;
      }}
      return out;
    }}
    L.r = de(L.r, L.rs);
    L.b = de(L.b, L.bs);
  }}

  var ALL_DETS = {all_dets};

  function unpackDets(D) {{
    if (D.rows) return;
    // dets is only present when the burst is missing a detector; the common
    // case is all fourteen in the standard order and is left implicit.
    if (!D.dets) D.dets = ALL_DETS;
    var u = b64bytes(D.d), nb = D.nb, rows = {{}};
    for (var i = 0; i < D.dets.length; i++) {{
      rows[D.dets[i]] = u.subarray(i * nb, (i + 1) * nb);
    }}
    D.rows = rows;
  }}

  // ---------------------------------------------------------- lightcurve
  function drawLC(r) {{
    var lcsvg = document.getElementById('lc');
    var L = r.L;
    if (!L) {{ lcsvg.innerHTML = ''; return; }}
    unpackLC(L);
    var COL = colors();
    var W = 960, H = 300, M = {{l: 60, r: 14, t: 12, b: 38}};
    var iw = W - M.l - M.r, ih = H - M.t - M.b;

    var e = L.e, rate = L.r, bkg = L.b;
    var x0 = e[0], x1 = e[e.length - 1];
    var lo = Math.min.apply(null, bkg), hi = Math.max.apply(null, rate);
    var head = Math.max((hi - lo) * 0.10, 1);
    var y0 = Math.max(0, lo - head), y1 = hi + head;

    function px(v) {{ return M.l + iw * (v - x0) / (x1 - x0 || 1); }}
    function py(v) {{ return M.t + ih * (y1 - v) / (y1 - y0 || 1); }}

    var out = [];

    function band(a, b, color, opacity, dashed) {{
      var xa = px(Math.max(a, x0)), xb = px(Math.min(b, x1));
      if (xb <= xa) return;
      // The outline is drawn at full strength rather than as a faint tint:
      // a low-alpha wash is legible on white and turns to mud on dark.
      out.push('<rect x="' + xa.toFixed(1) + '" y="' + M.t + '" width="' +
        (xb - xa).toFixed(1) + '" height="' + ih + '" fill="' + color +
        '" fill-opacity="' + opacity + '"' +
        (dashed ? ' stroke="' + color + '" stroke-width="1.4" stroke-dasharray="5,3"'
                : '') + ' />');
    }}

    // Background windows first, then the full T90 span outlined, then the
    // window actually used filled -- so what was discarded stays visible.
    (L.wbg || []).forEach(function (w) {{ band(w[0], w[1], COL.off, 0.16); }});
    band(L.w90[0], L.w90[1], COL.t90, 0.10, true);
    band(L.wsn[0], L.wsn[1], COL.t90, 0.34);

    out.push('<rect class="frame" x="' + M.l + '" y="' + M.t +
             '" width="' + iw + '" height="' + ih + '" />');

    function ticks(a, b, n) {{
      var raw = (b - a) / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
      var step = [1, 2, 2.5, 5, 10].map(function (m) {{ return m * mag; }})
                  .filter(function (v) {{ return v >= raw; }})[0] || mag * 10;
      var t = [], v = Math.ceil(a / step) * step;
      for (; v <= b; v += step) t.push(Math.round(v * 1e6) / 1e6);
      return t;
    }}

    out.push('<g class="axis">');
    ticks(x0, x1, 5).forEach(function (v) {{
      var X = px(v);
      out.push('<line x1="' + X + '" y1="' + M.t + '" x2="' + X + '" y2="' + (M.t + ih) + '"/>');
      out.push('<text x="' + X + '" y="' + (M.t + ih + 14) +
               '" font-size="9" text-anchor="middle">' + v + '</text>');
    }});
    ticks(y0, y1, 4).forEach(function (v) {{
      var Y = py(v);
      out.push('<line x1="' + M.l + '" y1="' + Y + '" x2="' + (M.l + iw) + '" y2="' + Y + '"/>');
      out.push('<text x="' + (M.l - 6) + '" y="' + (Y + 3) +
               '" font-size="9" text-anchor="end">' + v + '</text>');
    }});
    out.push('</g>');

    // Step plot: trigdat changes bin width across the burst, so the rate is
    // drawn across each bin rather than between bin centres.
    var d = '';
    for (var i = 0; i < rate.length; i++) {{
      var xa = px(e[i]), xb = px(e[i + 1]), Y = py(rate[i]);
      d += (i ? 'L' : 'M') + xa.toFixed(1) + ',' + Y.toFixed(1) +
           'L' + xb.toFixed(1) + ',' + Y.toFixed(1);
    }}
    out.push('<path d="' + d + '" fill="none" stroke="' + COL.off +
             '" stroke-width="1" stroke-opacity=".85" stroke-linejoin="round"/>');

    var db = '';
    for (var j = 0; j < bkg.length; j++) {{
      var xm = (px(e[j]) + px(e[j + 1])) / 2;
      db += (j ? 'L' : 'M') + xm.toFixed(1) + ',' + py(bkg[j]).toFixed(1);
    }}
    out.push('<path d="' + db + '" fill="none" stroke="' + COL.peak +
             '" stroke-width="1.6"/>');

    out.push('<text x="' + (M.l + iw / 2) + '" y="' + (H - 4) +
             '" font-size="10" text-anchor="middle">seconds from trigger</text>');
    out.push('<text transform="translate(12,' + (M.t + ih / 2) +
             ') rotate(-90)" font-size="10" text-anchor="middle">counts/s</text>');
    out.push('<text x="' + (M.l + 6) + '" y="' + (M.t + 12) +
             '" font-size="9">' + L.ts + ' ms bins, 50-300 keV</text>');
    lcsvg.innerHTML = out.join('');
  }}

  // --------------------------------------------------------------- panel
  function esc(s) {{
    return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) {{
      return {{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}}[c];
    }});
  }}

  function panel(r) {{
    var d = r.ps - r.as;
    var cls = r.ps < 2 ? 'good' : (r.ps > 10 ? 'bad' : '');
    var dcls = d < 0 ? 'good' : (d > 0 ? 'bad' : '');
    var h = '';
    // The panel spans the full page width, so readings are split into columns;
    // one long list would strand each label a thousand pixels from its value.
    h += '<div class="facthead">';
    h += '<h3>' + esc(r.t) + '</h3>';
    h += '<span class="sub">' + esc(r.n) + ' &middot; T90 ' + r.t90 + ' s (' + r.k +
         ') &middot; detectors ' + esc(r.det.join(', ')) + '</span>';
    h += '<span class="big ' + cls + '">' + r.ps.toFixed(2) + '\\u00B0</span>';
    h += '<span class="delta ' + dcls + '">' +
         (d === 0 ? 'unchanged from T90 window'
                  : (d < 0 ? 'improved by ' + Math.abs(d).toFixed(2) + '\\u00B0'
                           : 'worse by ' + d.toFixed(2) + '\\u00B0') + ' vs T90') +
         '</span>';
    h += '</div>';
    if (r.warn) {{
      h += '<p class="burstwarn"><strong>Flagged.</strong> ' + esc(r.warn) + '</p>';
    }}
    if (r.note) {{
      h += '<p class="burstnote"><strong>Note.</strong> ' + esc(r.note) + '</p>';
    }}
    h += '<div class="factgrid">';
    h += '<dl class="kv">'
      + '<dt>Official centroid</dt><dd>' + r.oc[0].toFixed(2) + ', ' + r.oc[1].toFixed(2) + '</dd>'
      + '<dt>Mine, T90 window</dt><dd>' + r.ac[0].toFixed(2) + ', ' + r.ac[1].toFixed(2) + '</dd>'
      + '<dt>Mine, best S/N</dt><dd>' + r.pc[0].toFixed(2) + ', ' + r.pc[1].toFixed(2) + '</dd>'
      + '<dt>Catalog position</dt><dd>' + r.cat[0].toFixed(2) + ', ' + r.cat[1].toFixed(2) + '</dd>'
      + '</dl>';
    h += '<dl class="kv">'
      + '<dt>Separation, T90</dt><dd>' + r.as.toFixed(2) + '\\u00B0</dd>'
      + '<dt>Separation, best S/N</dt><dd>' + r.ps.toFixed(2) + '\\u00B0</dd>'
      + '<dt>Window, T90</dt><dd>' + r.aw.toFixed(2) + ' s</dd>'
      + '<dt>Window, best S/N</dt><dd>' + r.pw.toFixed(2) + ' s</dd>'
      + '</dl>';
    h += '<dl class="kv">'
      + '<dt>90% area, official</dt><dd>' + r.oa + ' sq&deg;</dd>'
      + '<dt>90% area, mine</dt><dd>' + r.pa + ' sq&deg;</dd>'
      + '<dt>Binning</dt><dd>' + (r.L ? r.L.ts + ' ms' : '&mdash;') + '</dd>'
      + (r.ra_area ? '<dt>90% area, RoboBA kernel</dt><dd>' + r.ra_area +
         ' sq&deg;</dd>' : '')
      + '</dl>';
    h += '</div>';
    if (r.x) {{
      var best = Math.min(r.x.do, r.x.da, r.x.dp);
      function mark(v) {{ return v === best ? ' \\u2190 closest' : ''; }}
      h += '<div class="extbox wide">'
        + '<strong>' + esc(r.x.s) + '</strong> position ' + esc(r.x.n)
        + ' &mdash; ' + r.x.c[0].toFixed(4) + ', ' + r.x.c[1].toFixed(4)
        + '<dl class="kv" style="margin-top:.5rem">'
        + '<dt>Official is off by</dt><dd>' + r.x.do.toFixed(2) + '\\u00B0' + mark(r.x.do) + '</dd>'
        + '<dt>T90 window is off by</dt><dd>' + r.x.da.toFixed(2) + '\\u00B0' + mark(r.x.da) + '</dd>'
        + '<dt>Best S/N is off by</dt><dd>' + r.x.dp.toFixed(2) + '\\u00B0' + mark(r.x.dp) + '</dd>'
        + '</dl></div>';
    }} else {{
      h += '<p class="noext">No external position for this burst, so there is '
        + 'nothing here independent of GBM to check against.</p>';
    }}
    facts.innerHTML = h;
  }}

  // ------------------------------------------------------ detector panels
  // GBM has no optics. It works out where a burst was by comparing how much
  // each detector saw, so these fourteen panels are the actual input to the
  // fit rather than an illustration of it.
  function drawDets(r) {{
    var host = document.getElementById('dets');
    var D = r.D;
    if (!D) {{
      host.innerHTML = '<p class="panelcap">No TTE data for this burst.</p>';
      return;
    }}
    unpackDets(D);
    var COL = colors();
    var trg = {{}}, dead = {{}};
    (D.trg || []).forEach(function (d) {{ trg[d] = 1; }});
    (D.dead || []).forEach(function (d) {{ dead[d] = 1; }});

    var W = 150, H = 76, M = {{l: 3, r: 3, t: 13, b: 11}};
    var iw = W - M.l - M.r, ih = H - M.t - M.b;
    var span = D.hi - D.lo || 1;
    function px(t) {{ return M.l + iw * (t - D.lo) / span; }}

    function cell(det) {{
      var row = D.rows[det];
      if (!row) return '';
      var i = D.dets.indexOf(det), net = D.net[i];
      var bgo = det.charAt(0) === 'b';
      var sc = bgo ? D.sb : D.sn;
      var out = ['<svg class="dsvg" viewBox="0 0 ' + W + ' ' + H + '" role="img" ' +
                 'aria-label="' + det + ', net ' + net + ' counts">'];

      // Same window shading as the summed lightcurve above, so the eye can
      // carry across between the two.
      function band(a, b, color, op) {{
        var xa = px(Math.max(a, D.lo)), xb = px(Math.min(b, D.hi));
        if (xb <= xa) return;
        out.push('<rect x="' + xa.toFixed(1) + '" y="' + M.t + '" width="' +
                 (xb - xa).toFixed(1) + '" height="' + ih + '" fill="' + color +
                 '" fill-opacity="' + op + '" />');
      }}
      if (r.L) {{
        band(r.L.w90[0], r.L.w90[1], COL.t90, 0.13);
        band(r.L.wsn[0], r.L.wsn[1], COL.t90, 0.30);
      }}
      // The trigger instant, so each panel can be read without hunting for
      // the shared axis below it.
      if (D.lo < 0 && D.hi > 0) {{
        out.push('<line x1="' + px(0).toFixed(1) + '" y1="' + M.t + '" x2="' +
                 px(0).toFixed(1) + '" y2="' + (M.t + ih) + '" stroke="' +
                 COL.muted + '" stroke-width=".7" stroke-dasharray="2,2" ' +
                 'stroke-opacity=".8" />');
      }}
      out.push('<rect class="frame" x="' + M.l + '" y="' + M.t +
               '" width="' + iw + '" height="' + ih + '" />');

      var stroke = dead[det] ? COL.muted : (trg[det] ? COL.peak : COL.off);
      var d = '', w = iw / row.length;
      for (var k = 0; k < row.length; k++) {{
        var y = (M.t + ih * (1 - row[k] / 255)).toFixed(1);
        var xa = (M.l + k * w).toFixed(1), xb = (M.l + (k + 1) * w).toFixed(1);
        d += (k ? 'L' : 'M') + xa + ',' + y + 'L' + xb + ',' + y;
      }}
      out.push('<path d="' + d + '" fill="none" stroke="' + stroke +
               '" stroke-width="1" stroke-linejoin="round" />');
      out.push('<text x="' + M.l + '" y="9" font-size="9" class="dlab' +
               (trg[det] ? ' on' : '') + '">' + det +
               (trg[det] ? ' \\u2605' : '') + '</text>');
      out.push('<text x="' + (W - M.r) + '" y="9" font-size="8.5" ' +
               'text-anchor="end" class="dnet">' +
               (net > 0 ? '+' : '') + net.toLocaleString() + '</text>');
      out.push('<text x="' + (W - M.r) + '" y="' + (H - 2) + '" font-size="7.5" ' +
               'text-anchor="end" class="dsc">' + Math.round(sc[0]) + '-' +
               Math.round(sc[0] + sc[1]) + ' c/s</text>');
      out.push('</svg>');
      return out.join('');
    }}

    var nai = ['n0','n1','n2','n3','n4','n5','n6','n7','n8','n9','na','nb'];
    var html = '<div class="detgrid">' + nai.map(cell).join('') + '</div>';
    html += '<div class="detsplit">BGO &mdash; 200 keV to 40 MeV, and no part '
         +  'of the localization fit</div>';
    html += '<div class="detgrid bgo">' + ['b0','b1'].map(cell).join('') + '</div>';

    // These panels are zoomed on the burst, unlike the summed lightcurve
    // above which spans all of trigdat, so the range has to be stated.
    var bin = (D.hi - D.lo) / D.nb;
    var nd = bin < 0.1 ? 0 : (bin < 1 ? 1 : 2);
    html += '<p class="panelcap detaxis">' +
      D.lo.toFixed(nd < 2 ? 2 : 1) + ' to +' + D.hi.toFixed(nd < 2 ? 2 : 1) +
      ' s from trigger &middot; ' +
      (bin < 1 ? Math.round(bin * 1000) + ' ms' : bin.toFixed(2) + ' s') +
      ' bins &middot; dashed line is the trigger &middot; NaI {nai_band}, ' +
      'shown on a common scale</p>';
    host.innerHTML = html;
  }}

  // ------------------------------------------------- energy-resolved view
  var CH_BANDS = {ch_bands};
  var CH_EMIN = {ch_emin}, CH_EMAX = {ch_emax};

  // Perceptual ramp for the spectrogram, near-black at zero significance so
  // the quiet sky recedes. Interpolated through fixed stops rather than
  // computed from channel formulae: a hand-rolled version had blue rise then
  // drop at 30% of the scale, which put a false seam across the map and
  // broke the one property a colour scale has to have, that brighter means
  // more. Luminance here rises monotonically from 0.000 to 0.949.
  var RAMP = [[0,0,4],[16,8,47],[51,5,91],[87,16,110],[120,28,109],
              [154,41,98],[187,55,84],[217,74,63],[240,105,37],
              [251,148,16],[252,255,164]];

  function hot(f) {{
    f = (f < 0 ? 0 : (f > 1 ? 1 : f)) * (RAMP.length - 1);
    var i = Math.min(Math.floor(f), RAMP.length - 2), g = f - i;
    var a = RAMP[i], b = RAMP[i + 1];
    return 'rgb(' + Math.round(a[0] + (b[0] - a[0]) * g) + ',' +
                    Math.round(a[1] + (b[1] - a[1]) * g) + ',' +
                    Math.round(a[2] + (b[2] - a[2]) * g) + ')';
  }}

  // Band curves get a soft-to-hard blue-to-red progression instead of the
  // ramp above. No sub-range of a perceptual ramp reads well on both themes
  // as a thin line -- its bright end scores 1.05 against white -- whereas
  // these two endpoints are theme tokens already checked in both.
  function mix(c1, c2, f) {{
    function parse(h) {{
      h = (h || '').trim();
      return h.charAt(0) === '#' && h.length >= 7
        ? [parseInt(h.substr(1, 2), 16), parseInt(h.substr(3, 2), 16),
           parseInt(h.substr(5, 2), 16)]
        : null;
    }}
    var a = parse(c1), b = parse(c2);
    if (!a || !b) return c1;
    return 'rgb(' + Math.round(a[0] + (b[0] - a[0]) * f) + ',' +
                    Math.round(a[1] + (b[1] - a[1]) * f) + ',' +
                    Math.round(a[2] + (b[2] - a[2]) * f) + ')';
  }}

  function drawChan(r) {{
    var host = document.getElementById('chan');
    var C = r.C;
    if (!C) {{
      host.innerHTML = '<p class="panelcap">No TTE data for this burst.</p>';
      return;
    }}
    var COL = colors();
    var sg = b64bytes(C.sg), bd = b64bytes(C.bd);
    var nt = C.nt, nc = C.nc, nb = CH_BANDS.length - 1;

    // ---- spectrogram
    var W = 960, H = 300, M = {{l: 58, r: 92, t: 12, b: 34}};
    var iw = W - M.l - M.r, ih = H - M.t - M.b;
    var cw = iw / nt, chh = ih / nc;
    var out = ['<svg class="spec" viewBox="0 0 ' + W + ' ' + H + '" role="img" ' +
      'aria-label="Time against energy, shaded by significance above background">'];
    for (var c = 0; c < nc; c++) {{
      // Row 0 is the lowest energy, drawn at the bottom.
      var y = M.t + ih - (c + 1) * chh;
      for (var i = 0; i < nt; i++) {{
        var v = sg[c * nt + i];
        if (v < 6) continue;                     // leave the panel ground showing
        out.push('<rect x="' + (M.l + i * cw).toFixed(2) + '" y="' + y.toFixed(2) +
          '" width="' + (cw + 0.6).toFixed(2) + '" height="' + (chh + 0.6).toFixed(2) +
          '" fill="' + hot(v / 255) + '" />');
      }}
    }}
    out.push('<rect class="frame" x="' + M.l + '" y="' + M.t + '" width="' + iw +
             '" height="' + ih + '" />');

    function tx(t) {{ return M.l + iw * (t - C.lo) / (C.hi - C.lo || 1); }}
    if (r.L) {{
      [[r.L.w90, '5,3'], [r.L.wsn, null]].forEach(function (p) {{
        var a = tx(Math.max(p[0][0], C.lo)), b = tx(Math.min(p[0][1], C.hi));
        if (b <= a) return;
        out.push('<rect x="' + a.toFixed(1) + '" y="' + M.t + '" width="' +
          (b - a).toFixed(1) + '" height="' + ih + '" fill="none" stroke="' +
          COL.t90 + '" stroke-width="1.2"' +
          (p[1] ? ' stroke-dasharray="' + p[1] + '"' : '') + ' />');
      }});
    }}

    out.push('<g class="axis">');
    for (var c = 0; c <= nc; c += 4) {{
      var e = CH_EMIN * Math.pow(CH_EMAX / CH_EMIN, c / nc);
      var Y = M.t + ih - c * chh;
      out.push('<text x="' + (M.l - 7) + '" y="' + (Y + 3) + '" font-size="9" ' +
               'text-anchor="end">' + (e < 100 ? e.toFixed(0) : Math.round(e / 10) * 10) + '</text>');
    }}
    var span = C.hi - C.lo;
    for (var k = 0; k <= 4; k++) {{
      var tv = C.lo + span * k / 4, X = M.l + iw * k / 4;
      out.push('<text x="' + X.toFixed(1) + '" y="' + (M.t + ih + 14) +
               '" font-size="9" text-anchor="middle">' +
               (span > 40 ? tv.toFixed(0) : tv.toFixed(1)) + '</text>');
    }}
    out.push('<text x="' + (M.l + iw / 2) + '" y="' + (H - 4) +
             '" font-size="9" text-anchor="middle">seconds from trigger</text>');
    out.push('<text x="14" y="' + (M.t + ih / 2) + '" font-size="9" ' +
             'text-anchor="middle" transform="rotate(-90 14 ' +
             (M.t + ih / 2) + ')">keV</text>');
    out.push('</g>');

    // Colour key, as a strip of the same ramp.
    var kx = M.l + iw + 20, kw = 13, steps = 40;
    for (var s = 0; s < steps; s++) {{
      out.push('<rect x="' + kx + '" y="' + (M.t + ih - (s + 1) * ih / steps).toFixed(2) +
        '" width="' + kw + '" height="' + (ih / steps + 0.6).toFixed(2) +
        '" fill="' + hot(s / (steps - 1)) + '" />');
    }}
    out.push('<rect class="frame" x="' + kx + '" y="' + M.t + '" width="' + kw +
             '" height="' + ih + '" />');
    out.push('<g class="axis">');
    out.push('<text x="' + (kx + kw + 5) + '" y="' + (M.t + 8) + '" font-size="9">' +
             C.sr[1].toFixed(1) + '</text>');
    out.push('<text x="' + (kx + kw + 5) + '" y="' + (M.t + ih) + '" font-size="9">0</text>');
    out.push('<text x="' + (kx + kw + 5) + '" y="' + (M.t + ih / 2) +
             '" font-size="9">sigma</text>');
    out.push('</g></svg>');
    var html = out.join('');

    // ---- band lightcurves
    var BW = 150, BH = 76, BM = {{l: 3, r: 3, t: 13, b: 11}};
    var biw = BW - BM.l - BM.r, bih = BH - BM.t - BM.b;
    var cells = [];
    for (var j = 0; j < nb; j++) {{
      var g = ['<svg class="dsvg" viewBox="0 0 ' + BW + ' ' + BH + '" role="img" ' +
        'aria-label="' + CH_BANDS[j] + ' to ' + CH_BANDS[j + 1] + ' keV, net ' +
        C.tot[j] + ' counts">'];
      function bx(t) {{ return BM.l + biw * (t - C.lo) / (C.hi - C.lo || 1); }}
      if (r.L) {{
        [[r.L.w90, 0.13], [r.L.wsn, 0.30]].forEach(function (p) {{
          var a = bx(Math.max(p[0][0], C.lo)), b = bx(Math.min(p[0][1], C.hi));
          if (b > a) g.push('<rect x="' + a.toFixed(1) + '" y="' + BM.t + '" width="' +
            (b - a).toFixed(1) + '" height="' + bih + '" fill="' + COL.t90 +
            '" fill-opacity="' + p[1] + '" />');
        }});
      }}
      if (C.lo < 0 && C.hi > 0) {{
        g.push('<line x1="' + bx(0).toFixed(1) + '" y1="' + BM.t + '" x2="' +
          bx(0).toFixed(1) + '" y2="' + (BM.t + bih) + '" stroke="' + COL.muted +
          '" stroke-width=".7" stroke-dasharray="2,2" stroke-opacity=".8" />');
      }}
      g.push('<rect class="frame" x="' + BM.l + '" y="' + BM.t + '" width="' + biw +
             '" height="' + bih + '" />');
      var d = '', bw = biw / nt;
      for (var i = 0; i < nt; i++) {{
        var yy = (BM.t + bih * (1 - bd[j * nt + i] / 255)).toFixed(1);
        d += (i ? 'L' : 'M') + (BM.l + i * bw).toFixed(1) + ',' + yy +
             'L' + (BM.l + (i + 1) * bw).toFixed(1) + ',' + yy;
      }}
      g.push('<path d="' + d + '" fill="none" stroke="' +
             mix(COL.off, COL.peak, j / (nb - 1)) +
             '" stroke-width="1" stroke-linejoin="round" />');
      g.push('<text x="' + BM.l + '" y="9" font-size="8.5" class="dlab">' +
             CH_BANDS[j] + '-' + CH_BANDS[j + 1] + ' keV</text>');
      g.push('<text x="' + (BW - BM.r) + '" y="9" font-size="8.5" ' +
             'text-anchor="end" class="dnet">' + C.tot[j].toLocaleString() + '</text>');
      g.push('<text x="' + (BW - BM.r) + '" y="' + (BH - 2) + '" font-size="7.5" ' +
             'text-anchor="end" class="dsc">' + Math.round(C.bs[j][0]) + '-' +
             Math.round(C.bs[j][1]) + ' c/s</text>');
      g.push('</svg>');
      cells.push(g.join(''));
    }}
    html += '<div class="detsplit">The same counts split into eight bands. '
         +  'Each band is scaled to its own range, unlike the detector panels '
         +  'above &mdash; here the question is how the pulse SHAPE changes '
         +  'with energy, and a shared scale would flatten the hard bands to '
         +  'nothing, since 8-15 keV outcounts 600-900 keV many times over. '
         +  'Soft bands are drawn blue, hard ones red. Figures are net counts '
         +  'in T90, and the counts/s range each panel spans.</div>';
    html += '<div class="detgrid bands">' + cells.join('') + '</div>';
    html += '<p class="panelcap detaxis">summed over ' + C.dets.join(', ') +
            ' &middot; ' + Math.round((C.hi - C.lo) / nt * 1000) + ' ms bins ' +
            '&middot; 128 TTE channels rebinned to ' + nc + ' logarithmic ' +
            'channels, which is still finer than the detectors resolve</p>';
    host.innerHTML = html;
  }}

  function show() {{
    var rec = current().filter(function (r) {{ return r.t === grbSel.value; }})[0];
    if (!rec) return;
    draw(rec);
    drawAzim(rec);
    drawLC(rec);
    drawDets(rec);
    drawChan(rec);
    panel(rec);
  }}

  function step(n) {{
    var i = grbSel.selectedIndex + n;
    if (i >= 0 && i < grbSel.options.length) {{ grbSel.selectedIndex = i; show(); }}
  }}

  clsSel.addEventListener('change', function () {{ fillGrbs(null); show(); }});
  sortSel.addEventListener('change', function () {{ fillGrbs(grbSel.value); show(); }});
  grbSel.addEventListener('change', show);
  document.getElementById('prev').addEventListener('click', function () {{ step(-1); }});
  document.getElementById('next').addEventListener('click', function () {{ step(1); }});
  document.addEventListener('keydown', function (e) {{
    if (e.target.tagName === 'SELECT') return;
    if (e.key === 'ArrowLeft') step(-1);
    if (e.key === 'ArrowRight') step(1);
  }});

  fillGrbs(null);
  show();
}})();
</script>
"""

    with open(OUT, "w") as handle:
        handle.write(html)
    non_ascii = sum(1 for ch in html if ord(ch) > 127)
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1e6:.2f} MB, "
          f"{len(records)} bursts: {len(longs)} long / {len(shorts)} short, "
          f"non-ascii={non_ascii})")


if __name__ == "__main__":
    main()
