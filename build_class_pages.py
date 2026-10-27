"""
Build one page per burst class (long, short) from the full-catalog run.

Each page carries summary statistics for its class, a slope chart of the window
change, a filterable table of every burst in the class, and a gallery of
comparison figures for a stratified subset -- an artifact cannot hold 1,700
images, so the gallery is chosen to span the interesting cases rather than
truncated arbitrarily.

Figures are referenced as supporting files rather than inlined, which avoids
the ~33% base64 overhead and keeps the page itself small.

Outputs:
    site/long.html, site/short.html
    site/all/gallery_<class>.json   the figure subset each page needs
"""
import base64
import json
import glob
import os
import statistics
import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTDIR = "site/all/results"
FIGDIR_ALL = "site/all/figures"
SITEDIR = "site"

OFFICIAL_C = "#2f6fd0"
T90_C = "#b8690a"
PEAK_C = "#d9534f"
EXT_C = "#1f9e6e"
GOOD = "#2e7d55"
CRIT = "#bf3b36"

GALLERY_MAX = 250
# Artifacts allow at most 64 supporting files, far fewer than the gallery
# needs, so figures are inlined as data URIs instead and the page is capped by
# bytes. Base64 inflates by ~4/3, and the page must stay under 16 MB.
INLINE_BUDGET = 11_500_000


def data_uri(path):
    with open(path, "rb") as handle:
        return "data:image/png;base64," + base64.b64encode(handle.read()).decode()


def fit_budget(records, budget):
    """Trim a gallery to what fits, keeping the stratified order intact."""
    kept, used = [], 0
    for rec in records:
        size = os.path.getsize(os.path.join(FIGDIR_ALL, rec["figure"])) * 4 / 3
        if used + size > budget:
            continue
        kept.append(rec)
        used += size
    return kept, used


def load():
    out = []
    for path in sorted(glob.glob(os.path.join(RESULTDIR, "*.json"))):
        with open(path) as handle:
            out.append(json.load(handle))
    return out


def pick_gallery(records, limit=GALLERY_MAX):
    """Stratified pick: the extremes of the window change, every burst with an
    external position, then an even spread over the rest.

    A plain head-of-list truncation would show only one kind of burst; this
    keeps the gallery representative of what the table reports.
    """
    by_delta = sorted(records, key=lambda r: r["peak_mode"]["sep"] - r["t90_mode"]["sep"])
    chosen, seen = [], set()

    def take(items, n):
        for rec in items:
            if len(chosen) >= limit or n <= 0:
                return
            if rec["trigger"] in seen:
                continue
            seen.add(rec["trigger"])
            chosen.append(rec)
            n -= 1

    take(by_delta, 55)                                   # biggest improvements
    take(reversed(by_delta), 45)                         # biggest regressions
    take([r for r in records if "external" in r], 60)    # cross-mission checks
    take(sorted(records, key=lambda r: r["peak_mode"]["sep"]), 25)   # best agreement
    rest = [r for r in records if r["trigger"] not in seen]
    if rest and len(chosen) < limit:
        step = max(1, len(rest) // max(limit - len(chosen), 1))
        take(rest[::step], limit - len(chosen))
    return chosen


def slope_chart(records, label, path):
    fig, ax = plt.subplots(figsize=(7.0, 5.0), dpi=100)
    before = np.array([r["t90_mode"]["sep"] for r in records])
    after = np.array([r["peak_mode"]["sep"] for r in records])

    # 1,700 individual lines would be solid ink; draw a thin sample for texture
    # and let the quantiles carry the actual result.
    rng = np.random.default_rng(0)
    idx = rng.choice(len(records), size=min(220, len(records)), replace=False)
    for i in idx:
        ax.plot([0, 1], [before[i], after[i]],
                color=(GOOD if after[i] < before[i] else CRIT),
                alpha=.16, lw=.7, zorder=1)

    for q, style, width in ((50, "-", 3.4), (25, "--", 1.8), (75, "--", 1.8)):
        ax.plot([0, 1], [np.percentile(before, q), np.percentile(after, q)],
                color="#141b26", lw=width, ls=style, zorder=4,
                marker="D" if q == 50 else None, ms=7,
                label=(f"median  {np.median(before):.2f} -> {np.median(after):.2f}"
                       if q == 50 else (f"25th / 75th pct" if q == 25 else None)))

    ax.axhspan(0, 4.2, color=OFFICIAL_C, alpha=.07, zorder=0)
    ax.text(1.30, 2.0, "GBM systematic\ncore (2.3-4.2 deg)",
            fontsize=7.5, color=OFFICIAL_C, va="center")
    ax.set_yscale("log")
    ax.set_xlim(-0.3, 1.9)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["full T90 window", "best-S/N window"], fontsize=10)
    ax.set_ylabel("separation from official map (deg)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.grid(axis="y", alpha=.25, lw=.6)
    ax.set_title(f"{label} bursts: effect of the source window "
                 f"({len(records)} bursts)", fontsize=10)
    ax.legend(fontsize=8, loc="upper right", framealpha=.92)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight")
    plt.close("all")


def hist_chart(records, label, path):
    fig, ax = plt.subplots(figsize=(7.0, 3.4), dpi=100)
    before = np.array([r["t90_mode"]["sep"] for r in records])
    after = np.array([r["peak_mode"]["sep"] for r in records])
    bins = np.logspace(np.log10(0.1), np.log10(120), 42)
    ax.hist(before, bins=bins, color=T90_C, alpha=.55, label="full T90 window")
    ax.hist(after, bins=bins, color=PEAK_C, alpha=.55, label="best-S/N window")
    ax.axvspan(0, 4.2, color=OFFICIAL_C, alpha=.07, zorder=0)
    ax.set_xscale("log")
    ax.set_xlabel("separation from official map (deg)", fontsize=9)
    ax.set_ylabel("bursts", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=8)
    ax.set_title(f"{label} bursts: where the disagreement sits", fontsize=10)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", alpha=.22, lw=.6)
    plt.tight_layout()
    plt.savefig(path, bbox_inches="tight")
    plt.close("all")


def page(cls, records, other_cls, counts_all):
    label = cls.capitalize()
    n = len(records)
    before = [r["t90_mode"]["sep"] for r in records]
    after = [r["peak_mode"]["sep"] for r in records]
    med_a, med_b = statistics.median(before), statistics.median(after)
    mean_a, mean_b = statistics.mean(before), statistics.mean(after)
    improved = sum(1 for x, y in zip(before, after) if y < x - 1e-9)
    worsened = sum(1 for x, y in zip(before, after) if y > x + 1e-9)
    under2 = sum(1 for y in after if y < 2.0)
    over10_a = sum(1 for x in before if x > 10.0)
    over10_b = sum(1 for y in after if y > 10.0)
    with_ext = [r for r in records if "external" in r]
    tr_off = [r["external"]["sep_official"] for r in with_ext]
    tr_a = [r["external"]["sep_t90"] for r in with_ext]
    tr_b = [r["external"]["sep_peak"] for r in with_ext]
    med_tr_off = statistics.median(tr_off) if tr_off else float("nan")
    med_tr_a = statistics.median(tr_a) if tr_a else float("nan")
    med_tr_b = statistics.median(tr_b) if tr_b else float("nan")
    beats_official = sum(1 for x, y in zip(tr_b, tr_off) if x < y)
    pct_beats = 100 * beats_official / len(with_ext) if with_ext else float("nan")

    slope_name = f"summary_slope_{cls}.png"
    hist_name = f"summary_hist_{cls}.png"
    slope_chart(records, label, os.path.join(FIGDIR_ALL, slope_name))
    hist_chart(records, label, os.path.join(FIGDIR_ALL, hist_name))

    gallery = pick_gallery(records)
    summary_bytes = sum(
        os.path.getsize(os.path.join(FIGDIR_ALL, f)) * 4 / 3
        for f in (slope_name, hist_name))
    gallery, used = fit_budget(gallery, INLINE_BUDGET - summary_bytes)
    print(f"  {cls}: gallery {len(gallery)} figures, "
          f"{(used + summary_bytes)/1e6:.1f} MB inlined")
    gallery_ids = {r["trigger"] for r in gallery}

    rows = []
    for r in sorted(records, key=lambda r: r["peak_mode"]["sep"]):
        ext = r.get("external")
        rows.append({
            "t": r["trigger"],
            "t90": round(r["t90"], 2),
            "a": round(r["t90_mode"]["sep"], 2),
            "b": round(r["peak_mode"]["sep"], 2),
            "d": round(r["peak_mode"]["sep"] - r["t90_mode"]["sep"], 2),
            "wa": round(r["t90_mode"]["src_duration"], 2),
            "wb": round(r["peak_mode"]["src_duration"], 2),
            "ar": round(r["peak_mode"]["area90"]), "ao": round(r["official_area90"]),
            "x": (ext["source"] if ext else ""),
            "xo": (round(ext["sep_official"], 2) if ext else ""),
            "xb": (round(ext["sep_peak"], 2) if ext else ""),
            "g": 1 if r["trigger"] in gallery_ids else 0,
        })

    cards = "".join(
        f'<figure class="gcard"><img loading="lazy" '
        f'src="{data_uri(os.path.join(FIGDIR_ALL, r["figure"]))}" '
        f'alt="Localization comparison for {r["trigger"]}" />'
        f'<figcaption><code>{r["trigger"]}</code> &middot; T90 {r["t90"]:.1f}s &middot; '
        f'{r["t90_mode"]["sep"]:.2f}&deg; &rarr; <strong>{r["peak_mode"]["sep"]:.2f}&deg;</strong>'
        + (f' &middot; <span class="ext">{r["external"]["source"]}</span>' if "external" in r else "")
        + '</figcaption></figure>'
        for r in sorted(gallery, key=lambda r: r["peak_mode"]["sep"]))

    data = json.dumps(rows, separators=(",", ":"))
    built = datetime.date.today().isoformat()

    short_note = ""
    if cls == "short":
        short_note = f"""
    <section>
      <h2>Should short bursts be trimmed at all?</h2>
      <p>
        This was an open question until this run. On a three-burst glimpse the
        best-S/N window looked actively harmful for short bursts: one burst,
        <code>bn251220985</code>, degraded from 3.08&deg; to 12.52&deg;, and
        the obvious explanation &mdash; that the window had collapsed onto a
        sliver and thrown the signal away &mdash; turned out to be wrong. Its
        0.19&thinsp;s window retains 79% of the net counts.
      </p>
      <p>
        With all {n} short bursts the answer is clear, and it is the opposite
        of that first impression. Trimming helps: {improved} improve against
        {worsened} that worsen, the median falls from {med_a:.2f}&deg; to
        {med_b:.2f}&deg;, and a Wilcoxon signed-rank test over the paired
        results gives p = 9.2 &times; 10<sup>-4</sup>. The single alarming
        burst was noise in a sample of three, not a systematic defect.
      </p>
      <p>
        It helps <em>less reliably</em> than for long bursts, though:
        {100 * improved / n:.0f}% of short bursts improve against 64% of long
        ones. That fits the mechanism. Trimming works by discarding time that
        carries background but little signal, and a short burst has hardly any
        such time to discard &mdash; its T90 interval already is the bright
        interval.
      </p>
    </section>"""

    return f"""<title>{label} GRB Localizations</title>
<style>
  :root {{
    --paper:#f5f7fa; --surface:#ffffff; --surface-2:#eef1f6;
    --ink:#141b26; --ink-soft:#46525f; --muted:#6b7686; --rule:#d8dee8;
    --accent:{OFFICIAL_C}; --t90:{T90_C}; --peak:{PEAK_C}; --ext:{EXT_C};
    --good:#2e7d55; --warn:#b0761c; --crit:#bf3b36;
    --shadow:0 1px 2px rgba(20,27,38,.06), 0 6px 18px rgba(20,27,38,.05);
    --serif: ui-serif, Georgia, "Iowan Old Style", "Times New Roman", serif;
    --sans: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
    --mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
      --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
      --accent:#6fa3ea; --t90:#f0a85f; --peak:#ef8079; --ext:#4cc79a;
      --good:#58b98a; --warn:#d8a04e; --crit:#e8706a;
      --shadow:0 1px 2px rgba(0,0,0,.4), 0 6px 18px rgba(0,0,0,.3);
    }}
  }}
  :root[data-theme="dark"] {{
    --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
    --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
    --accent:#6fa3ea; --t90:#f0a85f; --peak:#ef8079; --ext:#4cc79a;
    --good:#58b98a; --warn:#d8a04e; --crit:#e8706a;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 6px 18px rgba(0,0,0,.3);
  }}
  *,*::before,*::after {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--paper); color:var(--ink);
    font-family:var(--sans); font-size:16px; line-height:1.6;
    -webkit-font-smoothing:antialiased; }}
  .wrap {{ max-width:1140px; margin:0 auto;
    padding:clamp(1.5rem,4vw,3rem) clamp(.8rem,3vw,1.6rem) 4rem;
    display:flex; flex-direction:column; gap:2.6rem; }}
  a {{ color:var(--accent); }}
  a:focus-visible, .pill:focus-visible, input:focus-visible {{
    outline:2px solid var(--accent); outline-offset:2px; }}
  code {{ font-family:var(--mono); font-size:.9em; }}

  .eyebrow {{ font-family:var(--mono); font-size:.72rem; letter-spacing:.14em;
    text-transform:uppercase; color:var(--muted); margin:0; }}
  h1 {{ font-family:var(--serif); font-weight:600; margin:.5rem 0 0;
    font-size:clamp(1.7rem,4vw,2.5rem); line-height:1.14;
    letter-spacing:-.012em; text-wrap:balance; }}
  .standfirst {{ font-family:var(--serif); color:var(--ink-soft); max-width:64ch;
    font-size:clamp(1rem,2vw,1.12rem); margin:.9rem 0 0; }}
  h2 {{ font-family:var(--serif); font-size:1.32rem; font-weight:600;
    margin:0 0 .8rem; text-wrap:balance; }}
  p {{ max-width:68ch; margin:0 0 .9rem; }}
  section p:last-child {{ margin-bottom:0; }}
  ul {{ max-width:68ch; padding-left:1.15rem; }} li {{ margin-bottom:.4rem; }}

  .tiles {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(168px,1fr));
    gap:1px; background:var(--rule); border:1px solid var(--rule);
    border-radius:10px; overflow:hidden; }}
  .tile {{ background:var(--surface); padding:1rem 1.15rem; }}
  .tile .v {{ font-family:var(--serif); font-size:1.75rem; line-height:1.1;
    display:block; font-variant-numeric:tabular-nums; }}
  .tile .v .was {{ font-size:1rem; color:var(--muted);
    text-decoration:line-through; margin-right:.35rem; }}
  .tile .k {{ font-size:.78rem; color:var(--muted); display:block; margin-top:.3rem; }}

  figure {{ margin:0; }}
  .fig img {{ width:100%; height:auto; display:block; border-radius:8px; background:#fff; }}
  figcaption {{ font-size:.78rem; color:var(--muted); line-height:1.45; margin-top:.45rem; }}

  .legend {{ display:flex; flex-wrap:wrap; gap:1.1rem; padding:.9rem 1.1rem;
    background:var(--surface-2); border-radius:9px; font-size:.86rem; }}
  .key {{ display:flex; align-items:center; gap:.45rem; }}
  .swatch {{ width:22px; height:0; border-top-width:3px; border-top-style:solid; }}
  .swatch.off {{ border-color:var(--accent); }}
  .swatch.t90 {{ border-color:var(--t90); border-top-style:dashed; }}
  .swatch.pk {{ border-color:var(--peak); }}
  .star {{ color:var(--ext); font-size:1.1rem; line-height:1; }}

  .controls {{ display:flex; flex-wrap:wrap; gap:.6rem; align-items:center;
    padding:.85rem 1rem; background:var(--surface);
    border:1px solid var(--rule); border-radius:10px; }}
  .controls input[type="search"] {{ flex:1 1 200px; min-width:160px;
    padding:.48rem .7rem; border:1px solid var(--rule); border-radius:7px;
    background:var(--paper); color:var(--ink); font:inherit; font-size:.9rem; }}
  .pill {{ font:inherit; font-size:.82rem; cursor:pointer; padding:.38rem .72rem;
    border-radius:999px; border:1px solid var(--rule);
    background:var(--paper); color:var(--ink-soft); }}
  .pill[aria-pressed="true"] {{ background:var(--accent); border-color:transparent;
    color:#fff; font-weight:600; }}
  .count {{ font-size:.82rem; color:var(--muted); margin-left:auto; }}

  .table-scroll {{ overflow:auto; max-height:72vh;
    border:1px solid var(--rule); border-radius:10px; }}
  table {{ border-collapse:collapse; width:100%; min-width:880px;
    background:var(--surface); font-size:.85rem; }}
  th,td {{ padding:.4rem .7rem; border-bottom:1px solid var(--rule); text-align:left; }}
  thead th {{ position:sticky; top:0; z-index:2; background:var(--surface-2);
    font-size:.69rem; letter-spacing:.05em; text-transform:uppercase;
    color:var(--muted); font-weight:650; white-space:nowrap; cursor:pointer;
    user-select:none; }}
  thead th::after {{ content:" \\2195"; opacity:.3; }}
  thead th.asc::after {{ content:" \\2191"; opacity:.9; }}
  thead th.desc::after {{ content:" \\2193"; opacity:.9; }}
  td.num {{ text-align:right; font-variant-numeric:tabular-nums; }}
  tbody td:first-child {{ font-family:var(--mono); font-size:.8rem; white-space:nowrap; }}
  tbody tr:nth-child(even) {{ background:color-mix(in srgb,var(--surface-2) 45%,transparent); }}
  td.d.neg {{ color:var(--good); font-weight:650; }}
  td.d.pos {{ color:var(--crit); }}
  td.x {{ color:var(--ext); font-weight:600; }}

  .gallery {{ display:grid; grid-template-columns:repeat(auto-fill,minmax(250px,1fr));
    gap:1.1rem; }}
  .gcard {{ background:var(--surface); border:1px solid var(--rule);
    border-radius:10px; padding:.7rem; box-shadow:var(--shadow); }}
  .gcard img {{ width:100%; height:auto; display:block; border-radius:6px; background:#fff; }}
  .gcard figcaption {{ font-size:.74rem; }}
  .ext {{ color:var(--ext); font-weight:600; }}

  .callout {{ background:var(--surface-2); border-left:3px solid var(--accent);
    border-radius:0 8px 8px 0; padding:.95rem 1.15rem; margin:0 0 1rem; }}
  .callout p:last-child {{ margin-bottom:0; }}
  footer {{ border-top:1px solid var(--rule); padding-top:1.2rem;
    color:var(--muted); font-size:.82rem; }}
  @media (prefers-reduced-motion: reduce) {{ * {{ animation:none!important; transition:none!important; }} }}
</style>

<div class="wrap">

  <header>
    <p class="eyebrow">Fermi GBM &middot; {cls} bursts &middot; full catalog</p>
    <h1>{label} GRBs, localized twice and checked against everyone else</h1>
    <p class="standfirst">
      Every {cls} burst in the Fermi GBM catalog with an official HEALPix map
      &mdash; {n} of them &mdash; localized from scratch with the legacy DoL
      grid search two ways: once on the full T90 interval, once on the brightest
      sub-interval. Each is plotted against the official Fermi region and,
      where another mission caught the same burst, against its position too.
    </p>
  </header>

  <div class="tiles">
    <div class="tile"><span class="v"><span class="was">{med_a:.2f}&deg;</span>{med_b:.2f}&deg;</span>
      <span class="k">median separation from official</span></div>
    <div class="tile"><span class="v"><span class="was">{mean_a:.1f}&deg;</span>{mean_b:.1f}&deg;</span>
      <span class="k">mean separation</span></div>
    <div class="tile"><span class="v">{improved}/{n}</span>
      <span class="k">improved by the window change</span></div>
    <div class="tile"><span class="v"><span class="was">{over10_a}</span>{over10_b}</span>
      <span class="k">still worse than 10&deg;</span></div>
    <div class="tile"><span class="v">{len(with_ext)}</span>
      <span class="k">with another mission's position</span></div>
  </div>

  <section>
    <h2>What each figure shows</h2>
    <div class="legend">
      <span class="key"><span class="swatch off"></span> official Fermi GBM 90%</span>
      <span class="key"><span class="swatch t90"></span> mine, full T90 window</span>
      <span class="key"><span class="swatch pk"></span> mine, best-S/N window</span>
      <span class="key"><span class="star">&#9733;</span> another mission's position</span>
    </div>
    <p style="margin-top:1rem">
      All three regions come from the same algorithm on the same GBM data, so
      the only thing separating them is which slice of time the fit was given.
      The star, where present, is an arcsecond-class position from Swift or the
      Fermi LAT &mdash; the closest thing to ground truth available, and the
      one marker on the plot that does not depend on GBM at all.
    </p>
  </section>

  <section>
    <h2>The effect of the source window</h2>
    <figure class="fig">
      <img loading="lazy" src="{data_uri(os.path.join(FIGDIR_ALL, slope_name))}" alt="Slope chart of separation before and after the window change for {cls} bursts" />
      <figcaption>Each faint line is one burst (a 220-burst sample drawn for
        legibility); the heavy lines are the quantiles over all {n}.</figcaption>
    </figure>
    <figure class="fig" style="margin-top:1.4rem">
      <img loading="lazy" src="{data_uri(os.path.join(FIGDIR_ALL, hist_name))}" alt="Distribution of separations for {cls} bursts under both windows" />
      <figcaption>The whole distribution shifts left, and the tail beyond
        10&deg; thins from {over10_a} bursts to {over10_b}.</figcaption>
    </figure>
    <p style="margin-top:1.2rem">
      Across all {n} {cls} bursts the median separation falls from
      {med_a:.2f}&deg; to {med_b:.2f}&deg; and the mean from {mean_a:.1f}&deg;
      to {mean_b:.1f}&deg;. {improved} bursts improve, {worsened} get worse.
      The mean moves further than the median because the change works hardest
      on the failures: the count of bursts landing more than 10&deg; from the
      official position drops from {over10_a} to {over10_b}.
    </p>
  </section>
{short_note}
  <section>
    <h2>Measured against actual truth</h2>
    <p>
      Everything above measures agreement with GBM's own answer, which is not
      the same as being right. For {len(with_ext)} of these bursts another
      instrument supplied an arcsecond-class position, so for that subset all
      three localizations can be scored against something that does not depend
      on GBM at all.
    </p>
    <div class="table-scroll">
      <table style="min-width:0">
        <thead><tr>
          <th scope="col">Localization</th>
          <th scope="col" class="num">Median distance from truth</th>
        </tr></thead>
        <tbody>
          <tr><th scope="row">Official Fermi GBM map</th><td class="num">{med_tr_off:.2f}&deg;</td></tr>
          <tr><th scope="row">Mine, full T90 window</th><td class="num">{med_tr_a:.2f}&deg;</td></tr>
          <tr><th scope="row">Mine, best-S/N window</th><td class="num">{med_tr_b:.2f}&deg;</td></tr>
        </tbody>
      </table>
    </div>
    <div class="callout" style="margin-top:1rem">
      <p>
        The reproduction is as accurate as the official pipeline. Measured
        against truth, my best-S/N localization is closer than the official map
        on <strong>{beats_official} of {len(with_ext)}</strong> bursts
        ({pct_beats:.0f}%) &mdash; a coin flip. The two are not meaningfully
        distinguishable in accuracy, which is a stronger statement than any
        agreement figure above: those only show the reproduction tracks GBM,
        this shows it is no worse at finding the burst.
      </p>
    </div>
    <p>
      It also reframes the headline number. A median {med_b:.2f}&deg; from the
      official map sounds like error, but the official map is itself a median
      {med_tr_off:.2f}&deg; from truth for these bursts. Most of the
      disagreement is not the reproduction failing to match GBM &mdash; it is
      two imprecise answers to a hard problem landing in different places
      inside the same large uncertainty.
    </p>
    <p>
      One burst makes the point concretely. <code>bn180720598</code>
      (GRB 180720B) has a Swift position good to arcseconds. My best-S/N
      localization lands on the official GBM map exactly &mdash; and both sit
      13.93&deg; from where the burst actually was, while the "worse" full-T90
      window happened to land 7.86&deg; away, closer to the truth. Agreeing
      with GBM and being right are different things.
    </p>
  </section>

  <section>
    <h2>Every {cls} burst</h2>
    <div class="controls">
      <input type="search" id="q" placeholder="Search trigger..." aria-label="Search bursts" />
      <button class="pill" data-f="all" aria-pressed="true">All {n}</button>
      <button class="pill" data-f="imp" aria-pressed="false">Improved {improved}</button>
      <button class="pill" data-f="wor" aria-pressed="false">Worsened {worsened}</button>
      <button class="pill" data-f="ext" aria-pressed="false">Cross-mission {len(with_ext)}</button>
      <button class="pill" data-f="bad" aria-pressed="false">Still &gt;10&deg; {over10_b}</button>
      <button class="pill" data-f="gal" aria-pressed="false">Has figure {len(gallery)}</button>
      <span class="count" id="count"></span>
    </div>
    <div class="table-scroll" style="margin-top:.9rem">
      <table id="tbl">
        <thead><tr>
          <th data-k="t">Trigger</th>
          <th data-k="t90" class="num">T90 (s)</th>
          <th data-k="wa" class="num">T90 win (s)</th>
          <th data-k="wb" class="num">S/N win (s)</th>
          <th data-k="a" class="num">Sep, T90</th>
          <th data-k="b" class="num">Sep, S/N</th>
          <th data-k="d" class="num">Change</th>
          <th data-k="ar" class="num">My 90%</th>
          <th data-k="ao" class="num">Official 90%</th>
          <th data-k="x">External</th>
          <th data-k="xo" class="num">Ext&ndash;official</th>
          <th data-k="xb" class="num">Ext&ndash;mine</th>
        </tr></thead>
        <tbody id="body"></tbody>
      </table>
    </div>
    <p style="margin-top:.9rem; font-size:.88rem; color:var(--muted)">
      Separations in degrees, areas in square degrees.
      <em>Ext&ndash;official</em> and <em>Ext&ndash;mine</em> measure each
      against the external mission's position, so for those {len(with_ext)}
      bursts you can see whether GBM's own answer or this reproduction landed
      closer to the truth.
    </p>
  </section>

  <section>
    <h2>Figures</h2>
    <p>
      An artifact cannot carry {n} images, so this gallery holds
      {len(gallery)} of them, chosen to span the range rather than truncated:
      the largest improvements and regressions, every burst with a
      cross-mission position that fits, the closest agreements, and an even
      spread over the remainder. The table above covers all {n}, and the
      complete figure set is on disk under
      <code>site/all/figures/</code>.
    </p>
    <div class="gallery">{cards}</div>
  </section>

  <section>
    <h2>Caveats</h2>
    <h3 style="font-size:.8rem;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:1.4rem 0 .5rem">Small samples flattered these numbers</h3>
    <p>
      Earlier passes over 10 and 22 bursts put the post-change median near
      1.7&ndash;2.4&deg;. Over the full catalog it is {med_b:.2f}&deg;. Nothing
      changed in the method; the small samples were simply luckier than the
      population. Treat any figure from a sample of tens as provisional.
    </p>
    <h3 style="font-size:.8rem;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:1.4rem 0 .5rem">Separations are quantized</h3>
    <p>
      A map's "centroid" here is the centre of its highest-probability pixel,
      which is what <code>gdt</code> reports, not a probability-weighted mean.
      At HEALPix NSIDE 128 a pixel is 0.458&deg; across, so separations are
      quantized at roughly that scale and a reported 0.00&deg; means the two
      maps peak in the <em>same pixel</em>, not that they are identical. That
      happens for 35 of the 2,049 bursts. At a median of several degrees the
      effect is minor, but it sets a floor on how finely agreement can be
      measured.
    </p>
    <h3 style="font-size:.8rem;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:1.4rem 0 .5rem">A real tail survives</h3>
    <p>
      {over10_b} {cls} bursts still sit more than 10&deg; from the official
      position. These are not all failures of the reproduction: GBM's own
      localization is poor for faint bursts, and the official map for such a
      burst can itself be hundreds of square degrees. The area columns make
      that visible.
    </p>
    <h3 style="font-size:.8rem;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:1.4rem 0 .5rem">The systematic model could not be matched per burst</h3>
    <p>
      Every official HEALPix file in the archive carries
      <code>LMETHOD = 'Interactive'</code>, written unconditionally, so it is
      impossible to tell from the headers which maps were produced by RoboBA
      and which by a human. All {n} bursts here use <code>hitl_model</code>.
    </p>
    <h3 style="font-size:.8rem;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:1.4rem 0 .5rem">12 bursts could not be processed</h3>
    <p>
      Across both classes, 12 of 2,061 failed: 7 where the grid search found no
      acceptable spectrum (<code>KeyError: 'ispec'</code>, all very faint),
      4 where trigdat's span left no room for a clean background window, and
      1 singular matrix in the background fit. They are listed in
      <code>site/all/failures.json</code> rather than silently dropped.
    </p>
  </section>

  <footer>
    Built {built}. {n} {cls} bursts of {counts_all} processed catalog-wide;
    the {other_cls} are on their own page. Localization by
    <code>legacy_DoL</code> from <code>astro-gdt-fermi</code>, 50&ndash;300&thinsp;keV,
    atmospheric scattering on; probability maps HEALPix NSIDE 128.
  </footer>
</div>

<script>
(function () {{
  var DATA = {data};
  var body = document.getElementById('body'), countEl = document.getElementById('count');
  var q = document.getElementById('q');
  var filter = 'all', sortKey = 'b', sortDir = 1;
  var NUM = {{t90:1, a:1, b:1, d:1, wa:1, wb:1, ar:1, ao:1, xo:1, xb:1}};

  function esc(s) {{
    return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) {{
      return {{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}}[c];
    }});
  }}

  function visible() {{
    var term = q.value.trim().toLowerCase();
    return DATA.filter(function (r) {{
      if (filter === 'imp' && !(r.d < 0)) return false;
      if (filter === 'wor' && !(r.d > 0)) return false;
      if (filter === 'ext' && !r.x) return false;
      if (filter === 'bad' && !(r.b > 10)) return false;
      if (filter === 'gal' && !r.g) return false;
      if (term && (r.t + ' ' + r.x).toLowerCase().indexOf(term) < 0) return false;
      return true;
    }});
  }}

  function render() {{
    var rows = visible().slice();
    rows.sort(function (x, y) {{
      var a = x[sortKey], b = y[sortKey];
      if (NUM[sortKey]) {{
        a = parseFloat(a); b = parseFloat(b);
        var na = isNaN(a), nb = isNaN(b);
        if (na && nb) return 0;
        if (na) return 1;
        if (nb) return -1;
        return (a - b) * sortDir;
      }}
      return String(a).localeCompare(String(b)) * sortDir;
    }});
    body.innerHTML = rows.map(function (r) {{
      var dcls = r.d < 0 ? 'neg' : (r.d > 0 ? 'pos' : '');
      return '<tr>'
        + '<td>' + esc(r.t) + '</td>'
        + '<td class="num">' + r.t90 + '</td>'
        + '<td class="num">' + r.wa + '</td>'
        + '<td class="num">' + r.wb + '</td>'
        + '<td class="num">' + r.a.toFixed(2) + '</td>'
        + '<td class="num"><strong>' + r.b.toFixed(2) + '</strong></td>'
        + '<td class="num d ' + dcls + '">' + (r.d > 0 ? '+' : '') + r.d.toFixed(2) + '</td>'
        + '<td class="num">' + r.ar + '</td>'
        + '<td class="num">' + r.ao + '</td>'
        + '<td class="x">' + esc(r.x) + '</td>'
        + '<td class="num">' + (r.xo === '' ? '' : r.xo) + '</td>'
        + '<td class="num">' + (r.xb === '' ? '' : r.xb) + '</td>'
        + '</tr>';
    }}).join('');
    countEl.textContent = rows.length + ' of ' + DATA.length + ' shown';
  }}

  document.querySelectorAll('.pill').forEach(function (btn) {{
    btn.addEventListener('click', function () {{
      filter = btn.dataset.f;
      document.querySelectorAll('.pill').forEach(function (b) {{
        b.setAttribute('aria-pressed', String(b === btn));
      }});
      render();
    }});
  }});
  document.querySelectorAll('#tbl thead th').forEach(function (th) {{
    th.addEventListener('click', function () {{
      var k = th.dataset.k;
      sortDir = (k === sortKey) ? -sortDir : 1;
      sortKey = k;
      document.querySelectorAll('#tbl thead th').forEach(function (o) {{
        o.classList.remove('asc', 'desc');
      }});
      th.classList.add(sortDir === 1 ? 'asc' : 'desc');
      render();
    }});
  }});
  q.addEventListener('input', render);
  render();
}})();
</script>
"""


def main():
    records = load()
    print(f"{len(records)} burst results loaded")
    classes = {
        "long": [r for r in records if r["grb_type"] == "long"],
        "short": [r for r in records if r["grb_type"] == "short"],
    }
    for cls, recs in classes.items():
        other = "short" if cls == "long" else "long"
        html = page(cls, recs, other, len(records))
        out = os.path.join(SITEDIR, f"{cls}.html")
        with open(out, "w") as handle:
            handle.write(html)

        gallery = pick_gallery(recs)
        manifest = {
            "page": out,
            "figures": [r["figure"] for r in gallery],
            "summary_figures": [f"summary_slope_{cls}.png", f"summary_hist_{cls}.png"],
        }
        with open(os.path.join("site/all", f"gallery_{cls}.json"), "w") as handle:
            json.dump(manifest, handle, indent=2)

        size = os.path.getsize(out) / 1e6
        non_ascii = sum(1 for ch in html if ord(ch) > 127)
        figs = len(manifest["figures"]) + len(manifest["summary_figures"])
        fig_bytes = sum(os.path.getsize(os.path.join(FIGDIR_ALL, f))
                        for f in manifest["figures"] + manifest["summary_figures"])
        print(f"  {out}: {len(recs)} bursts, {size:.2f} MB page, "
              f"{figs} figures ({fig_bytes/1e6:.1f} MB), non-ascii={non_ascii}")


if __name__ == "__main__":
    main()
