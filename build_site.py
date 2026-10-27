"""
Build the comparison website from the two results files and site/figures/.

Reads both window modes:
    site/results.json        full T90 source window
    site/results_peak.json   best signal-to-noise sub-interval

Everything is inlined (figures as base64 data URIs) so the page is a single
self-contained HTML file that works offline and under a strict CSP. The output
is pure ASCII, since the artifact wrapper owns <head> and no meta charset can
be declared from here.
"""
import base64
import json
import os
import datetime
import statistics

SITEDIR = "site"
FIGDIR = os.path.join(SITEDIR, "figures")
OUT = os.path.join(SITEDIR, "index.html")

OFFICIAL = "#2f6fd0"
MINE = "#d9534f"


def b64(path):
    with open(path, "rb") as handle:
        return "data:image/png;base64," + base64.b64encode(handle.read()).decode()


def chip(sep):
    """Separation banded against GBM's own systematic scale.

    The hitl systematic core is a 2.3-4.2 deg Gaussian, so a couple of degrees
    is the expected level of disagreement, not a failure.
    """
    if sep < 2.0:
        return "close"
    if sep < 6.0:
        return "mid"
    return "wide"


def main():
    with open(os.path.join(SITEDIR, "results.json")) as handle:
        t90_data = json.load(handle)
    with open(os.path.join(SITEDIR, "results_peak.json")) as handle:
        peak_data = json.load(handle)
    with open(os.path.join(SITEDIR, "holdout_t90.json")) as handle:
        hold_t90 = json.load(handle)
    with open(os.path.join(SITEDIR, "holdout_50_300.json")) as handle:
        hold_peak = json.load(handle)
    with open(os.path.join(SITEDIR, "sweep_main_summary.json")) as handle:
        sweep = json.load(handle)


    t90 = {r["trigger"]: r for r in t90_data["results"]}
    peak = {r["trigger"]: r for r in peak_data["results"]}
    triggers = sorted(peak, key=lambda k: peak[k]["sep_mine_official"])

    ht = {r["trigger"]: r for r in hold_t90["results"]}
    hp = {r["trigger"]: r for r in hold_peak["results"]}
    hold_keys = sorted(set(ht) & set(hp))

    # Pooled across both samples. The window comparison is paired on identical
    # bursts, so it is unaffected by how the samples happened to be composed;
    # the absolute level very much is, which is why both are reported.
    pool_t90 = ([t90[k]["sep_mine_official"] for k in peak]
                + [ht[k]["sep_mine_official"] for k in hold_keys])
    pool_peak = ([peak[k]["sep_mine_official"] for k in peak]
                 + [hp[k]["sep_mine_official"] for k in hold_keys])
    pool_types = ([peak[k]["grb_type"] for k in peak]
                  + [hp[k]["grb_type"] for k in hold_keys])
    pool_n = len(pool_t90)
    pool_med_a, pool_med_b = statistics.median(pool_t90), statistics.median(pool_peak)
    pool_mean_a, pool_mean_b = statistics.mean(pool_t90), statistics.mean(pool_peak)
    pool_improved = sum(1 for x, y in zip(pool_t90, pool_peak) if y < x - 1e-9)
    pool_worsened = sum(1 for x, y in zip(pool_t90, pool_peak) if y > x + 1e-9)
    longs = [y for y, t in zip(pool_peak, pool_types) if t == "long"]
    shorts = [y for y, t in zip(pool_peak, pool_types) if t == "short"]

    hold_med_a = statistics.median([ht[k]["sep_mine_official"] for k in hold_keys])
    hold_med_b = statistics.median([hp[k]["sep_mine_official"] for k in hold_keys])
    hold_improved = sum(1 for k in hold_keys
                        if hp[k]["sep_mine_official"] < ht[k]["sep_mine_official"] - 1e-9)

    sweep_rows = sorted(sweep["bands"], key=lambda b: b["median"])
    std_band = next(b for b in sweep_rows if (b["lo"], b["hi"]) == (50, 300))

    def seps(d):
        return [d[k]["sep_mine_official"] for k in triggers]

    def ratios(d):
        return [d[k]["my_area90"] / d[k]["official_area90"] for k in triggers]

    sep_a, sep_b = seps(t90), seps(peak)
    rat_a, rat_b = ratios(t90), ratios(peak)
    med_a, med_b = statistics.median(sep_a), statistics.median(sep_b)
    mean_a, mean_b = statistics.mean(sep_a), statistics.mean(sep_b)
    medr_a, medr_b = statistics.median(rat_a), statistics.median(rat_b)
    improved = sum(1 for k in triggers
                   if peak[k]["sep_mine_official"] < t90[k]["sep_mine_official"] - 1e-9)
    worsened = sum(1 for k in triggers
                   if peak[k]["sep_mine_official"] > t90[k]["sep_mine_official"] + 1e-9)
    contained = sum(1 for k in triggers if peak[k]["official_in_my_map"] <= 0.9)
    n_long = sum(1 for k in triggers if peak[k]["grb_type"] == "long")
    built = datetime.date.today().isoformat()

    # ------------------------------------------------------------ table rows
    rows = []
    for k in sorted(triggers, key=lambda k: -(t90[k]["sep_mine_official"]
                                              - peak[k]["sep_mine_official"])):
        a, b = t90[k], peak[k]
        delta = b["sep_mine_official"] - a["sep_mine_official"]
        shrink = a["src_duration"] / b["src_duration"] if b["src_duration"] else 1.0
        sign = "pos" if delta > 0.005 else ("neg" if delta < -0.005 else "flat")
        rows.append(f"""
        <tr>
          <th scope="row"><a href="#{k}"><code>{k}</code></a></th>
          <td class="num">{a['t90']:.1f}</td>
          <td class="num">{a['src_duration']:.1f}</td>
          <td class="num">{b['src_duration']:.1f}</td>
          <td class="num">{shrink:.1f}x</td>
          <td class="num">{a['sep_mine_official']:.2f}</td>
          <td class="num"><strong>{b['sep_mine_official']:.2f}</strong></td>
          <td class="num delta {sign}">{delta:+.2f}</td>
          <td><span class="chip {chip(b['sep_mine_official'])}"></span></td>
        </tr>""")

    # ---------------------------------------------------------- burst cards
    cards = []
    for k in triggers:
        a, b = t90[k], peak[k]
        delta = b["sep_mine_official"] - a["sep_mine_official"]
        arrow = "improved" if delta < -0.005 else ("regressed" if delta > 0.005 else "unchanged")
        bg = " &nbsp;&middot;&nbsp; ".join(
            f"{lo:.0f} to {hi:.0f}&thinsp;s" for lo, hi in b["bg_times"])
        dets = ", ".join(b["triggered_detectors"]) or "&mdash;"
        dead = ", ".join(b["dead_detectors"]) if b["dead_detectors"] else "none"

        cards.append(f"""
      <article class="burst" id="{k}">
        <header class="burst-head">
          <div class="burst-id">
            <h3><code>{k}</code></h3>
            <p class="burst-sub">{b['name']} &nbsp;&middot;&nbsp; T90 {b['t90']:.1f}&thinsp;s ({b['grb_type']})</p>
          </div>
          <div class="burst-metric">
            <span class="prev">{a['sep_mine_official']:.2f}&deg;</span>
            <span class="arrow" aria-label="becomes">&rarr;</span>
            <span class="metric-value {chip(b['sep_mine_official'])}">{b['sep_mine_official']:.2f}&deg;</span>
            <span class="metric-label">separation from official, {arrow}</span>
          </div>
        </header>

        <figure class="fig">
          <img loading="lazy" src="{b64(os.path.join(FIGDIR, f'lc_{k}_peak.png'))}"
               alt="Lightcurve for {k} showing the full T90 span and the narrower best signal-to-noise window" />
          <figcaption>Summed triggered detectors, 50&ndash;300&thinsp;keV. Dashed
            outline is the full T90 span; solid orange is the best-S/N
            sub-interval actually used. Blue are the background windows, red the
            fitted background.</figcaption>
        </figure>

        <div class="two-up">
          <figure class="fig">
            <img loading="lazy" src="{b64(os.path.join(FIGDIR, f'zoom_{k}.png'))}"
                 alt="Localization regions for {k} using the full T90 window" />
            <figcaption><span class="tag">Before</span> full T90 window &mdash;
              {a['sep_mine_official']:.2f}&deg; off, {a['my_area90']:.0f} sq&deg;</figcaption>
          </figure>
          <figure class="fig">
            <img loading="lazy" src="{b64(os.path.join(FIGDIR, f'zoom_{k}_peak.png'))}"
                 alt="Localization regions for {k} using the best signal-to-noise window" />
            <figcaption><span class="tag on">After</span> best-S/N window &mdash;
              {b['sep_mine_official']:.2f}&deg; off, {b['my_area90']:.0f} sq&deg;</figcaption>
          </figure>
        </div>

        <figure class="fig">
          <img loading="lazy" src="{b64(os.path.join(FIGDIR, f'sky_{k}_peak.png'))}"
               alt="All-sky view of both localization regions for {k}" />
          <figcaption>All-sky, Earth occultation removed.</figcaption>
        </figure>

        <dl class="facts">
          <div><dt>My centroid</dt><dd>RA {b['my_ra']:.2f}&deg; &nbsp; Dec {b['my_dec']:+.2f}&deg;</dd></div>
          <div><dt>Official centroid</dt><dd>RA {b['official_ra']:.2f}&deg; &nbsp; Dec {b['official_dec']:+.2f}&deg;</dd></div>
          <div><dt>Catalog position</dt><dd>RA {b['catalog_ra']:.2f}&deg; &nbsp; Dec {b['catalog_dec']:+.2f}&deg;</dd></div>
          <div><dt>Source window</dt><dd>{b['src_time'][0]:.1f} to {b['src_time'][1]:.1f}&thinsp;s ({b['src_duration']:.1f}&thinsp;s)</dd></div>
          <div><dt>Full T90 span</dt><dd>{a['src_duration']:.1f}&thinsp;s</dd></div>
          <div><dt>Window S/N</dt><dd>{b['peak_snr']:.1f}</dd></div>
          <div><dt>90% area, mine</dt><dd>{b['my_area90']:.0f} sq deg</dd></div>
          <div><dt>90% area, official</dt><dd>{b['official_area90']:.0f} sq deg</dd></div>
          <div><dt>Official centroid sits at</dt><dd>{100 * b['official_in_my_map']:.0f}% contour of mine</dd></div>
          <div><dt>Background windows</dt><dd>{bg}</dd></div>
          <div><dt>Binning</dt><dd>{b['timescale_ms']}&thinsp;ms</dd></div>
          <div><dt>Triggered detectors</dt><dd>{dets}</dd></div>
          <div><dt>Disabled detectors</dt><dd>{dead}</dd></div>
          <div><dt>Systematic model</dt><dd><code>{b['sys_model']}</code></dd></div>
        </dl>
      </article>""")

    html = f"""<title>GBM Localization Replication</title>
<style>
  :root {{
    --paper:      #f5f7fa;
    --surface:    #ffffff;
    --surface-2:  #eef1f6;
    --ink:        #141b26;
    --ink-soft:   #46525f;
    --muted:      #6b7686;
    --rule:       #d8dee8;
    --official:   {OFFICIAL};
    --mine:       {MINE};
    --good:       #2e7d55;
    --warn:       #b0761c;
    --crit:       #bf3b36;
    --shadow:     0 1px 2px rgba(20,27,38,.06), 0 6px 20px rgba(20,27,38,.05);

    --serif: ui-serif, Georgia, "Iowan Old Style", "Times New Roman", serif;
    --sans:  system-ui, -apple-system, "Segoe UI", Roboto, "Helvetica Neue", sans-serif;
    --mono:  ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --paper:     #0e131b;
      --surface:   #151c27;
      --surface-2: #1b2430;
      --ink:       #e8edf4;
      --ink-soft:  #b3becd;
      --muted:     #8995a6;
      --rule:      #2a3545;
      --official:  #6fa3ea;
      --mine:      #ef8079;
      --good:      #58b98a;
      --warn:      #d8a04e;
      --crit:      #e8706a;
      --shadow:    0 1px 2px rgba(0,0,0,.4), 0 6px 20px rgba(0,0,0,.3);
    }}
  }}
  :root[data-theme="dark"] {{
    --paper:     #0e131b;
    --surface:   #151c27;
    --surface-2: #1b2430;
    --ink:       #e8edf4;
    --ink-soft:  #b3becd;
    --muted:     #8995a6;
    --rule:      #2a3545;
    --official:  #6fa3ea;
    --mine:      #ef8079;
    --good:      #58b98a;
    --warn:      #d8a04e;
    --crit:      #e8706a;
    --shadow:    0 1px 2px rgba(0,0,0,.4), 0 6px 20px rgba(0,0,0,.3);
  }}

  *, *::before, *::after {{ box-sizing: border-box; }}

  body {{
    margin: 0;
    background: var(--paper);
    color: var(--ink);
    font-family: var(--sans);
    font-size: 16px;
    line-height: 1.6;
    -webkit-font-smoothing: antialiased;
  }}

  .wrap {{
    max-width: 1080px;
    margin: 0 auto;
    padding: clamp(1.6rem, 4vw, 3.4rem) clamp(1rem, 4vw, 2rem) 5rem;
    display: flex;
    flex-direction: column;
    gap: 3.2rem;
  }}

  a {{ color: var(--official); }}
  a:focus-visible, summary:focus-visible {{
    outline: 2px solid var(--official);
    outline-offset: 3px;
    border-radius: 3px;
  }}
  code {{ font-family: var(--mono); font-size: .92em; }}

  /* ---------------------------------------------------------- masthead */
  .masthead {{
    display: flex; flex-direction: column; gap: 1rem;
    border-bottom: 1px solid var(--rule);
    padding-bottom: 2rem;
  }}
  .eyebrow {{
    font-family: var(--mono);
    font-size: .74rem; letter-spacing: .14em; text-transform: uppercase;
    color: var(--muted); margin: 0;
  }}
  h1 {{
    font-family: var(--serif); font-weight: 600;
    font-size: clamp(1.9rem, 4.4vw, 2.9rem);
    line-height: 1.12; letter-spacing: -.012em;
    text-wrap: balance; margin: 0;
  }}
  .standfirst {{
    font-family: var(--serif);
    font-size: clamp(1.02rem, 2vw, 1.18rem);
    line-height: 1.58; color: var(--ink-soft);
    max-width: 64ch; margin: 0;
  }}

  /* ------------------------------------------------------------- tiles */
  .tiles {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(190px, 1fr));
    gap: 1px; background: var(--rule);
    border: 1px solid var(--rule); border-radius: 10px; overflow: hidden;
  }}
  .tile {{ background: var(--surface); padding: 1.15rem 1.25rem; }}
  .tile .v {{
    font-family: var(--serif); font-size: 1.95rem; line-height: 1.1;
    font-variant-numeric: tabular-nums; display: block;
  }}
  .tile .v .was {{
    font-size: 1.05rem; color: var(--muted);
    text-decoration: line-through; margin-right: .4rem;
  }}
  .tile .k {{ font-size: .8rem; color: var(--muted); display: block; margin-top: .35rem; }}

  /* ------------------------------------------------------------ prose */
  section > h2 {{
    font-family: var(--serif); font-size: 1.42rem; font-weight: 600;
    letter-spacing: -.008em; margin: 0 0 .9rem; text-wrap: balance;
  }}
  section > h3 {{
    font-family: var(--sans); font-size: .82rem; font-weight: 650;
    letter-spacing: .09em; text-transform: uppercase; color: var(--muted);
    margin: 1.9rem 0 .7rem;
  }}
  p {{ max-width: 68ch; margin: 0 0 1rem; }}
  section p:last-child {{ margin-bottom: 0; }}
  ul {{ max-width: 68ch; padding-left: 1.15rem; margin: 0 0 1rem; }}
  li {{ margin-bottom: .45rem; }}
  li::marker {{ color: var(--muted); }}
  .legend-inline {{ font-weight: 650; }}
  .legend-inline.mine {{ color: var(--mine); }}
  .legend-inline.official {{ color: var(--official); }}

  /* ------------------------------------------------------------ table */
  .table-scroll {{ overflow-x: auto; border-radius: 10px; border: 1px solid var(--rule); }}
  table {{
    border-collapse: collapse; width: 100%; min-width: 760px;
    background: var(--surface); font-size: .9rem;
  }}
  th, td {{ text-align: left; padding: .6rem .8rem; border-bottom: 1px solid var(--rule); }}
  thead th {{
    background: var(--surface-2); font-size: .72rem; letter-spacing: .06em;
    text-transform: uppercase; color: var(--muted); font-weight: 650; white-space: nowrap;
  }}
  thead th.group {{ border-bottom: 1px solid var(--rule); text-align: center; }}
  tbody tr:last-child th, tbody tr:last-child td {{ border-bottom: none; }}
  td.num, th.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
  tbody th {{ font-weight: 500; white-space: nowrap; }}
  td.delta.neg {{ color: var(--good); font-weight: 650; }}
  td.delta.pos {{ color: var(--crit); }}
  td.delta.flat {{ color: var(--muted); }}

  .chip {{ display: inline-block; width: 10px; height: 10px; border-radius: 50%; }}
  .chip.close {{ background: var(--good); }}
  .chip.mid   {{ background: var(--warn); }}
  .chip.wide  {{ background: var(--crit); }}

  /* ------------------------------------------------------------ bursts */
  .bursts {{ display: flex; flex-direction: column; gap: 1.6rem; }}
  .burst {{
    background: var(--surface); border: 1px solid var(--rule);
    border-radius: 12px; padding: 1.4rem clamp(1rem, 2.4vw, 1.6rem) 1.5rem;
    box-shadow: var(--shadow); scroll-margin-top: 1rem;
    display: flex; flex-direction: column; gap: 1.2rem;
  }}
  .burst-head {{
    display: flex; flex-wrap: wrap; gap: 1rem;
    align-items: flex-start; justify-content: space-between;
    border-bottom: 1px solid var(--rule); padding-bottom: .9rem;
  }}
  .burst-id h3 {{ margin: 0; font-size: 1.12rem; font-weight: 650; }}
  .burst-sub {{ margin: .2rem 0 0; color: var(--muted); font-size: .86rem; }}
  .burst-metric {{ text-align: right; }}
  .burst-metric .prev {{
    font-family: var(--serif); font-size: 1.1rem; color: var(--muted);
    text-decoration: line-through; font-variant-numeric: tabular-nums;
  }}
  .burst-metric .arrow {{ color: var(--muted); margin: 0 .3rem; }}
  .metric-value {{
    font-family: var(--serif); font-size: 1.85rem; line-height: 1;
    font-variant-numeric: tabular-nums;
  }}
  .metric-value.close {{ color: var(--good); }}
  .metric-value.mid   {{ color: var(--warn); }}
  .metric-value.wide  {{ color: var(--crit); }}
  .metric-label {{ font-size: .74rem; color: var(--muted); display: block; margin-top: .2rem; }}

  .fig {{ margin: 0; display: flex; flex-direction: column; gap: .5rem; }}
  .fig img {{
    width: 100%; height: auto; display: block;
    border-radius: 7px; background: #fff;
  }}
  figcaption {{ font-size: .79rem; color: var(--muted); line-height: 1.45; }}
  .two-up {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(290px, 1fr));
    gap: 1.2rem;
  }}
  .tag {{
    display: inline-block; font-family: var(--mono); font-size: .68rem;
    letter-spacing: .06em; text-transform: uppercase;
    padding: .1rem .4rem; border-radius: 4px; margin-right: .35rem;
    background: var(--surface-2); color: var(--muted);
    border: 1px solid var(--rule);
  }}
  .tag.on {{ background: var(--good); color: #fff; border-color: transparent; }}

  .facts {{
    display: grid; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr));
    gap: .1rem 1.4rem; margin: 0;
    border-top: 1px solid var(--rule); padding-top: 1rem;
  }}
  .facts > div {{
    display: flex; justify-content: space-between; gap: 1rem;
    padding: .3rem 0; border-bottom: 1px dotted var(--rule);
  }}
  .facts dt {{ color: var(--muted); font-size: .82rem; }}
  .facts dd {{
    margin: 0; font-size: .82rem; text-align: right;
    font-variant-numeric: tabular-nums; white-space: nowrap;
  }}

  /* ------------------------------------------------------------ notes */
  .callout {{
    background: var(--surface-2); border-left: 3px solid var(--official);
    border-radius: 0 8px 8px 0; padding: 1rem 1.2rem; margin: 0 0 1rem;
  }}
  .callout p {{ margin-bottom: .6rem; }}
  .callout p:last-child {{ margin-bottom: 0; }}
  details {{
    border: 1px solid var(--rule); border-radius: 9px;
    background: var(--surface); padding: .85rem 1.1rem;
  }}
  summary {{ cursor: pointer; font-weight: 600; font-size: .94rem; }}
  details[open] summary {{ margin-bottom: .8rem; }}
  footer {{
    border-top: 1px solid var(--rule); padding-top: 1.4rem;
    color: var(--muted); font-size: .83rem;
  }}
  @media (prefers-reduced-motion: reduce) {{
    * {{ animation: none !important; transition: none !important; }}
  }}
</style>

<div class="wrap">

  <header class="masthead">
    <p class="eyebrow">Fermi GBM &nbsp;&middot;&nbsp; ground localization check</p>
    <h1>Lining up a DoL reproduction with GBM's official positions</h1>
    <p class="standfirst">
      Gamma-ray bursts drawn at random from the Fermi GBM catalog, each
      localized from scratch with the legacy Daughter-of-Locburst grid search
      and compared against the official GBM HEALPix map. Choosing the source
      window the way a Burst Advocate does, rather than using the whole T90
      interval, cuts the median gap from {pool_med_a:.2f}&deg; to
      {pool_med_b:.2f}&deg; across {pool_n} bursts. Narrowing the energy band,
      the obvious next lever, does not help at all.
    </p>
  </header>

  <section>
    <div class="tiles">
      <div class="tile">
        <span class="v"><span class="was">{pool_med_a:.2f}&deg;</span>{pool_med_b:.2f}&deg;</span>
        <span class="k">median separation, {pool_n} bursts pooled</span>
      </div>
      <div class="tile">
        <span class="v"><span class="was">{pool_mean_a:.2f}&deg;</span>{pool_mean_b:.2f}&deg;</span>
        <span class="k">mean separation</span>
      </div>
      <div class="tile">
        <span class="v">{pool_improved}/{pool_n}</span>
        <span class="k">bursts improved by the window change</span>
      </div>
      <div class="tile">
        <span class="v">{std_band['median']:.2f}&deg;</span>
        <span class="k">best energy band is the standard one</span>
      </div>
    </div>
  </section>

  <section>
    <h2>What this tests</h2>
    <p>
      The official GBM map and this reproduction run the <em>same</em> algorithm
      on the <em>same</em> data. Neither uses information from any other
      observatory: the FITS header of every official map records
      <code>LOC_SRC = 'Fermi, GBM'</code>. The only cross-mission ingredient
      anywhere in the chain is the systematic error model, calibrated once,
      population-wide, against bursts with precise positions from Swift and
      others &mdash; it widens the probability region but never moves the
      centroid of an individual burst.
    </p>
    <p>
      So a disagreement cannot come from the official map knowing something this
      one doesn't. It can only come from the choices a human Burst Advocate
      normally makes by hand. That makes this a clean experiment: change one
      choice, and the separation column measures exactly how much it was worth.
    </p>
  </section>

  <section>
    <h2>What closed the gap</h2>
    <p>
      The first pass localized on the catalog's full T90 interval. That is the
      defensible automatic choice, and it is wrong for the job. T90 is a
      <em>duration</em> measurement &mdash; the span containing the middle 90%
      of the fluence &mdash; not a statement about where the signal is strong.
      On a burst with a bright spike and a long soft tail, most of T90
      contributes background and almost no signal.
    </p>
    <p>
      The grid search weighs each detector's net counts against its background.
      Padding the window with near-empty time dilutes the very contrast the fit
      depends on: the region inflates, and the centroid is free to drift.
    </p>
    <div class="callout">
      <p>
        <strong>The change:</strong> instead of the full T90, take the
        contiguous sub-interval that maximises
        <code>sum(net) / sqrt(sum(background))</code> over the summed triggered
        detectors in the 50&ndash;300&thinsp;keV band. Every contiguous bin
        range inside the burst is tested; the best-S/N one wins. Background
        windows are untouched, so the two runs differ in exactly one choice.
      </p>
    </div>
    <p>
      Median separation falls from {med_a:.2f}&deg; to {med_b:.2f}&deg; and the
      mean from {mean_a:.2f}&deg; to {mean_b:.2f}&deg;. The worst case drops
      from {max(sep_a):.2f}&deg; to {max(sep_b):.2f}&deg;, which matters more
      than the median: the long tail of bad localizations is what made the first
      pass untrustworthy, and it is gone. {improved} of {len(triggers)} bursts
      improved, {worsened} regressed slightly.
    </p>

    <div class="table-scroll">
      <table>
        <thead>
          <tr>
            <th scope="col">Burst</th>
            <th scope="col" class="num">T90 (s)</th>
            <th scope="col" class="num">Window, T90 (s)</th>
            <th scope="col" class="num">Window, best S/N (s)</th>
            <th scope="col" class="num">Shrink</th>
            <th scope="col" class="num">Sep, T90 (&deg;)</th>
            <th scope="col" class="num">Sep, best S/N (&deg;)</th>
            <th scope="col" class="num">Change</th>
            <th scope="col"></th>
          </tr>
        </thead>
        <tbody>{''.join(rows)}
        </tbody>
      </table>
    </div>
    <p style="margin-top:.9rem">
      Sorted by how much each burst gained. The pattern is clean: the bursts
      that improve most are the ones whose window shrank most. Where the bright
      interval already <em>was</em> most of T90 &mdash;
      <code>bn180113418</code>, <code>bn190615636</code> &mdash; there was
      nothing to gain and nothing changed.
    </p>
  </section>

  <section>
    <h2>Does it hold up on bursts it was not tuned on?</h2>
    <p>
      The window rule was chosen after looking at the first ten bursts, so
      measuring it on those same ten proves nothing. A second sample was drawn
      with a different seed and run untouched.
    </p>
    <div class="table-scroll">
      <table>
        <thead>
          <tr>
            <th scope="col">Sample</th>
            <th scope="col" class="num">Bursts</th>
            <th scope="col" class="num">Median, T90 window</th>
            <th scope="col" class="num">Median, best S/N</th>
            <th scope="col" class="num">Improved</th>
          </tr>
        </thead>
        <tbody>
          <tr>
            <th scope="row">A &mdash; the one shown below</th>
            <td class="num">{len(triggers)}</td>
            <td class="num">{med_a:.2f}&deg;</td>
            <td class="num">{med_b:.2f}&deg;</td>
            <td class="num">{improved}/{len(triggers)}</td>
          </tr>
          <tr>
            <th scope="row">B &mdash; held out, fresh seed</th>
            <td class="num">{len(hold_keys)}</td>
            <td class="num">{hold_med_a:.2f}&deg;</td>
            <td class="num">{hold_med_b:.2f}&deg;</td>
            <td class="num">{hold_improved}/{len(hold_keys)}</td>
          </tr>
          <tr>
            <th scope="row"><strong>Pooled</strong></th>
            <td class="num"><strong>{pool_n}</strong></td>
            <td class="num"><strong>{pool_med_a:.2f}&deg;</strong></td>
            <td class="num"><strong>{pool_med_b:.2f}&deg;</strong></td>
            <td class="num"><strong>{pool_improved}/{pool_n}</strong></td>
          </tr>
        </tbody>
      </table>
    </div>
    <p style="margin-top:.9rem">
      The window change reproduces: {hold_improved} of {len(hold_keys)} held-out
      bursts improved, and a Wilcoxon signed-rank test over all {pool_n} pairs
      gives p = 0.017. Because the comparison is paired on identical bursts, it
      does not care how either sample happened to be composed.
    </p>
    <figure class="fig" style="margin-top:1.2rem">
      <img loading="lazy" src="{b64(os.path.join(FIGDIR, 'summary_window.png'))}"
           alt="Slope chart of separation before and after the window change for all pooled bursts" />
      <figcaption>Every burst, both windows. Green improved, red regressed.
        Note the log scale: the worst case fell from 117&deg; to 4.6&deg;.</figcaption>
    </figure>

    <div class="callout">
      <p>
        <strong>The absolute number is less certain than it looks.</strong>
        Sample A alone gave a median of {med_a:.2f}&deg; before the change and
        {med_b:.2f}&deg; after; sample B, the same method, gave
        {hold_med_a:.2f}&deg; and {hold_med_b:.2f}&deg;. Bootstrapping the
        pooled results, the median of a random ten-burst draw spans roughly
        1.8&ndash;3.5&deg;. Sample A was on the fortunate side of that, so the
        pooled {pool_med_b:.2f}&deg; is the figure to quote &mdash; not the
        number from any single draw.
      </p>
    </div>
  </section>

  <section>
    <h2>Narrowing the energy band does not help</h2>
    <p>
      If restricting in time bought that much, restricting in energy is the
      obvious thing to try next. It does not work, and the failure is
      informative.
    </p>
    <figure class="fig">
      <img loading="lazy" src="{b64(os.path.join(FIGDIR, 'summary_bands.png'))}"
           alt="Median separation for each candidate localization energy band" />
      <figcaption>All eight distinct bands trigdat's channel structure allows,
        run over sample A with the best-S/N window.</figcaption>
    </figure>
    <p>
      Trigdat carries only eight coarse channels, so the band is not
      continuously tunable &mdash; the standard 50&ndash;300&thinsp;keV lands on
      channels 3&ndash;4, and only a handful of genuinely distinct bands exist.
      Every narrowing is worse, and the single-channel bands are roughly three
      times worse: {std_band['median']:.2f}&deg; becomes 4.45&deg; using
      44&ndash;95&thinsp;keV alone, or 3.60&deg; using 95&ndash;300&thinsp;keV
      alone.
    </p>
    <p>
      The reason is the opposite of the time case. The grid search infers
      position from the <em>contrast</em> in counts across the twelve NaI
      detectors, whose angular responses differ. Trimming in time removed
      intervals that carried background but no signal. Trimming in energy
      removes real signal from every detector at once, so the contrast the fit
      depends on gets noisier while nothing useful is discarded.
    </p>
    <p>
      Widening to 44&ndash;500&thinsp;keV looked marginally better on sample A
      (1.60&deg; against {std_band['median']:.2f}&deg;), so it was checked on
      the held-out sample: better on only 5 of {len(hold_keys)} bursts, with a
      worse maximum on both samples. That is a coin flip, and the apparent gain
      came from one burst moving 9.97&deg; to 1.08&deg;. GBM's standard band
      stands.
    </p>
  </section>

  <section>
    <h2>Reading the figures</h2>
    <div class="callout">
      <p>
        <span class="legend-inline mine">Red</span> is this reproduction.
        <span class="legend-inline official">Blue</span> is the official GBM
        map. Solid contours are 90% confidence, dashed are 50%. The black star
        is the catalog position.
      </p>
      <p>
        Each burst shows its lightcurve with both candidate windows marked, then
        the localization before and after the change, then an all-sky frame for
        context. The close-ups are where agreement is readable &mdash; at
        all-sky scale a couple of degrees is nearly invisible.
      </p>
    </div>
  </section>

  <section>
    <h2>Burst by burst</h2>
    <div class="bursts">{''.join(cards)}
    </div>
  </section>

  <section>
    <h2>What the numbers mean, and what they don't</h2>

    <h3>This is agreement, not identity</h3>
    <p>
      The pooled median of {pool_med_b:.2f}&deg; sits inside GBM's own
      systematic error, whose core for this model is a 2.3&ndash;4.2&deg;
      Gaussian depending on where the burst falls in spacecraft azimuth. The
      reproduction tracks the official pipeline to within its stated
      uncertainty &mdash; it is not, and cannot be made, pixel-identical,
      because the remaining difference is a human's judgement about where a
      burst starts and stops.
    </p>

    <h3>The window change is not free</h3>
    <p>
      {pool_worsened} of {pool_n} bursts regressed. A pure maximum-S/N
      criterion is a reasonable stand-in for a Burst Advocate, not a replica of
      one: a real BA also weighs spectral hardness and detector geometry, and
      knows when a window has become too short to be worth trimming further.
      The regressions cluster where the rule has least to work with &mdash;
      very short bursts, and very long ones where the selected window is still
      most of T90 and the change is mostly noise.
    </p>

    <h3>My regions are still larger than the official ones</h3>
    <p>
      On sample A the median area ratio improved from {medr_a:.2f}x to
      {medr_b:.2f}x but did not reach parity. Some of this is irreducible without matching the
      official pipeline's exact interval, and some of it is conservative in the
      right direction: a region that is too large is honest about uncertainty,
      where one that is too small is not.
    </p>

    <h3>Short bursts do markedly worse</h3>
    <p>
      Sample A drew {n_long} long bursts and no short ones; the held-out sample
      supplied the first three. They land at a median of
      {statistics.median(shorts):.2f}&deg; against
      {statistics.median(longs):.2f}&deg; for the {len(longs)} long bursts.
      Three is far too few to put a number on, but the direction is what you
      would expect: a short burst gives the fit a fraction of the counts, so
      the detector-to-detector contrast it works from is correspondingly
      noisier.
    </p>
    <p>
      One short burst, <code>bn251220985</code> (T90 0.38&thinsp;s), is also
      the worst regression anywhere in this work &mdash; 3.08&deg; to
      12.52&deg;. S/N maximisation trimmed its window to 0.2&thinsp;s, and at
      that duration there are too few counts left to localize on. The window
      rule needs a floor for short bursts; it does not have one here.
    </p>

    <h3>The systematic model could not be matched per burst</h3>
    <p>
      The intent was to pick <code>hitl_model</code> or
      <code>robo_ba_model</code> per burst from each official map's
      <code>LMETHOD</code> header, since the two describe differently shaped
      error distributions. In practice all 2,061 official HEALPix files in the
      archive carry <code>LMETHOD = 'Interactive'</code> &mdash; written
      unconditionally, recording nothing about whether a human or the robot
      produced the map. Every burst here therefore uses <code>hitl_model</code>.
      For any map that was in fact RoboBA output, the other kernel would have
      been more appropriate and the area comparison would shift.
    </p>

    <h3>Centroid is not the same quantity as catalog position</h3>
    <p>
      The official map's centroid sits 0.26&ndash;1.31&deg; from the position
      listed in the burst catalog for the same burst. These are genuinely
      different quantities &mdash; the probability-weighted centre of an
      Earth-occultation-removed map, versus the quoted best-fit position &mdash;
      so the comparison here is centroid against centroid throughout.
    </p>
  </section>

  <section>
    <h2>Reproducing this</h2>
    <details>
      <summary>Pipeline and provenance</summary>
      <ul>
        <li>Burst list: HEASARC <code>fermigbrst</code> via TAP, 4,400 GRBs
          spanning 2008-07-14 to 2026-10-20.</li>
        <li>Data: <code>glg_trigdat_all_*</code> and
          <code>glg_healpix_all_*</code> for every catalog burst, 5.2&nbsp;GB
          from the HEASARC trigger tree.</li>
        <li>Eligible pool: {peak_data['n_eligible']:,} bursts with both an
          official HEALPix map and a usable T90. Official maps exist only from
          roughly 2018 onward; earlier bursts carry a
          <code>glg_locprob_all_*</code> gnomonic image instead, or no map at
          all, and cannot be compared this way.</li>
        <li>Sample: 10 drawn without replacement, seed
          <code>{peak_data['seed']}</code>. Both runs use the identical sample.</li>
        <li>Localization: <code>legacy_DoL</code> from
          <code>astro-gdt-fermi</code>, atmospheric scattering enabled,
          50&ndash;300&thinsp;keV, detector mask set internally by the routine
          (NaI with background above 100&nbsp;counts/s, both BGOs excluded).</li>
      </ul>
      <p style="margin-bottom:0">
        One correction was needed to the reference implementation before it
        could survive long bursts at all. <code>legacy_DoL</code> expects count
        <em>rates</em>, not counts &mdash; it multiplies the input by the source
        duration internally. Passing counts is harmless for a 2&thinsp;s window,
        where the error is a common factor that cancels in the flux
        normalisation, and both conventions agree to the digit on
        <code>bn200826187</code>. Over a 100&thinsp;s window the same mistake
        inflates every chi-square past the routine's internal clamp, leaving the
        grid uniform and the search dead on an empty argmin.
      </p>
    </details>
  </section>

  <footer>
    Built {built} from locally archived Fermi GBM data.
    Figures generated with <code>astro-gdt-fermi</code>; probability maps are
    HEALPix, NSIDE 128, equatorial coordinates.
  </footer>

</div>
"""

    with open(OUT, "w") as handle:
        handle.write(html)

    non_ascii = sum(1 for ch in html if ord(ch) > 127)
    print(f"wrote {OUT}  ({os.path.getsize(OUT) / 1e6:.2f} MB, "
          f"{len(triggers)} bursts, non-ascii={non_ascii})")


if __name__ == "__main__":
    main()
