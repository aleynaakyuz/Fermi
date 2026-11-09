"""
Render the SAA background measurement as a standalone page.

Separate from the burst viewer on purpose. The viewer is trigger-centred and
already sits at 15.4 MB against a 16 MB cap; this is continuous-time data
answering a different question, and it fits in about a tenth of a megabyte.

Input:  site/saa_profiles.json   (extract_saa.py)
Output: site/saa.html
"""
import datetime
import json
import os

SRC = "site/saa_profiles.json"
OUT = "site/saa.html"

ACCENT = "#2f6fd0"
HOT = "#d9534f"
WARM = "#b8690a"
COOL = "#1f9e6e"

NAI = [f"n{d}" for d in "0123456789ab"]
BGO = ["b0", "b1"]


def path_from(xs, ys, px, py):
    """Polyline path, skipping gaps where a value is missing."""
    out, pen = [], False
    for x, y in zip(xs, ys):
        if y is None:
            pen = False
            continue
        out.append(("L" if pen else "M") + f"{px(x):.1f},{py(y):.1f}")
        pen = True
    return "".join(out)


def place_labels(entries, top, bottom, gap=14.0):
    """Nudge right-hand series labels apart so they stay readable.

    The twelve NaI curves have all but converged by 90 minutes, so their
    labels land within a pixel or two of each other and overlap. Sort by the
    position each wants, push each down to clear the one above, then shift
    the whole stack up if it has run past the bottom.

    entries: list of (wanted_y, payload). Returns [(placed_y, payload)].
    """
    if not entries:
        return []
    order = sorted(entries, key=lambda e: e[0])
    placed, y = [], None
    for want, payload in order:
        y = want if y is None else max(want, y + gap)
        placed.append([y, payload])
    overflow = placed[-1][0] - bottom
    if overflow > 0:
        for item in placed:
            item[0] -= overflow
    # And if that pushed the top one above the frame, compress downward again.
    if placed[0][0] < top:
        shift = top - placed[0][0]
        for item in placed:
            item[0] = min(item[0] + shift, bottom)
    return [(y, payload) for y, payload in placed]


def axis_ticks(lo, hi, n=5):
    import math
    raw = (hi - lo) / n
    if raw <= 0:
        return [lo]
    mag = 10 ** math.floor(math.log10(raw))
    step = next((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw), mag * 10)
    out, v = [], math.ceil(lo / step) * step
    while v <= hi + 1e-9:
        out.append(round(v, 6))
        v += step
    return out


def chart_day(day):
    """24 hours of one detector, with SAA passages marked."""
    W, H, M = 980, 260, dict(l=62, r=16, t=14, b=38)
    iw, ih = W - M["l"] - M["r"], H - M["t"] - M["b"]
    xs, ys = day["t"], day["r"]
    x0, x1 = 0.0, max(xs)
    # ys carries a None at each SAA break, so the pen lifts; skip those when
    # working out the range.
    finite = [v for v in ys if v is not None]
    lo, hi = min(finite), max(finite)
    pad = (hi - lo) * 0.08 or 1
    y0, y1 = max(0, lo - pad), hi + pad
    px = lambda v: M["l"] + iw * (v - x0) / (x1 - x0 or 1)
    py = lambda v: M["t"] + ih * (y1 - v) / (y1 - y0 or 1)

    s = [f'<svg class="fig" viewBox="0 0 {W} {H}" role="img" aria-label="'
         f'One day of count rate with South Atlantic Anomaly passages marked">']
    for a, b in day["saa"]:
        s.append(f'<rect class="saaband" x="{px(a):.1f}" y="{M["t"]}" '
                 f'width="{max(px(b)-px(a),1):.1f}" height="{ih}" />')
    s.append(f'<rect class="frame" x="{M["l"]}" y="{M["t"]}" width="{iw}" '
             f'height="{ih}" />')
    s.append('<g class="axis">')
    for v in axis_ticks(x0, x1, 8):
        if not (x0 <= v <= x1):
            continue
        s.append(f'<line x1="{px(v):.1f}" y1="{M["t"]}" x2="{px(v):.1f}" '
                 f'y2="{M["t"]+ih}"/>')
        s.append(f'<text x="{px(v):.1f}" y="{M["t"]+ih+14}" text-anchor="middle">'
                 f'{v:g}</text>')
    for v in axis_ticks(y0, y1, 4):
        if not (y0 <= v <= y1):
            continue
        s.append(f'<line x1="{M["l"]}" y1="{py(v):.1f}" x2="{M["l"]+iw}" '
                 f'y2="{py(v):.1f}"/>')
        s.append(f'<text x="{M["l"]-7}" y="{py(v)+3:.1f}" text-anchor="end">'
                 f'{v:g}</text>')
    s.append(f'<text x="{M["l"]+iw/2}" y="{H-4}" text-anchor="middle">'
             f'hours into {day["day"]} UTC</text>')
    s.append(f'<text x="13" y="{M["t"]+ih/2}" text-anchor="middle" '
             f'transform="rotate(-90 13 {M["t"]+ih/2})">counts/s</text>')
    s.append('</g>')
    s.append(f'<path class="line" d="{path_from(xs, ys, px, py)}" />')
    s.append('</svg>')
    return "".join(s)


def chart_fold(data, dets, title_colour):
    """Folded decay, several detectors overlaid."""
    W, H, M = 980, 320, dict(l=62, r=120, t=14, b=38)
    iw, ih = W - M["l"] - M["r"], H - M["t"] - M["b"]
    t = data["t"]
    series = [(d, data["fold"][d]) for d in dets if d in data["fold"]]
    vals = [v for _, ys in series for v in ys if v is not None]
    lo, hi = min(vals), max(vals)
    pad = (hi - lo) * 0.08 or 1
    y0, y1 = lo - pad, hi + pad
    px = lambda v: M["l"] + iw * v / (max(t) or 1)
    py = lambda v: M["t"] + ih * (y1 - v) / (y1 - y0 or 1)

    s = [f'<svg class="fig" viewBox="0 0 {W} {H}" role="img" aria-label="'
         f'Count rate against time since leaving the SAA">']
    s.append(f'<rect class="frame" x="{M["l"]}" y="{M["t"]}" width="{iw}" '
             f'height="{ih}" />')
    s.append('<g class="axis">')
    for v in axis_ticks(0, max(t), 6):
        if v > max(t):
            continue
        s.append(f'<line x1="{px(v):.1f}" y1="{M["t"]}" x2="{px(v):.1f}" '
                 f'y2="{M["t"]+ih}"/>')
        s.append(f'<text x="{px(v):.1f}" y="{M["t"]+ih+14}" text-anchor="middle">'
                 f'{v:g}</text>')
    for v in axis_ticks(y0, y1, 4):
        if not (y0 <= v <= y1):
            continue
        s.append(f'<line x1="{M["l"]}" y1="{py(v):.1f}" x2="{M["l"]+iw}" '
                 f'y2="{py(v):.1f}"/>')
        s.append(f'<text x="{M["l"]-7}" y="{py(v)+3:.1f}" text-anchor="end">'
                 f'{v:g}</text>')
    s.append(f'<text x="{M["l"]+iw/2}" y="{H-4}" text-anchor="middle">'
             f'minutes since leaving the SAA</text>')
    s.append(f'<text x="13" y="{M["t"]+ih/2}" text-anchor="middle" '
             f'transform="rotate(-90 13 {M["t"]+ih/2})">counts/s</text>')
    s.append('</g>')
    wanted = []
    for i, (det, ys) in enumerate(series):
        frac = i / max(len(series) - 1, 1)
        opacity = 0.35 + 0.6 * (1 - frac)
        s.append(f'<path class="line" style="stroke:{title_colour};'
                 f'stroke-opacity:{opacity:.2f}" '
                 f'd="{path_from(t, ys, px, py)}" />')
        last = next((v for v in reversed(ys) if v is not None), None)
        if last is not None:
            wanted.append((py(last), (det, title_colour, 0.45 + 0.55 * (1 - frac))))
    for y, (det, colour, op) in place_labels(wanted, M["t"] + 11, M["t"] + ih - 2):
        s.append(f'<text class="serieslab" x="{M["l"]+iw+6}" y="{y+3:.1f}" '
                 f'style="fill:{colour};fill-opacity:{op:.2f}">{det}</text>')
    s.append('</svg>')
    return "".join(s)


def chart_chan(chan, t):
    """Energy dependence: each channel relative to its own quiet level."""
    W, H, M = 980, 320, dict(l=62, r=150, t=14, b=38)
    iw, ih = W - M["l"] - M["r"], H - M["t"] - M["b"]
    rows = chan["rows"]
    vals = [v for r in rows for v in r["rel"] if v is not None]
    y0, y1 = min(vals) - 0.03, max(vals) + 0.03
    px = lambda v: M["l"] + iw * v / (max(t) or 1)
    py = lambda v: M["t"] + ih * (y1 - v) / (y1 - y0 or 1)

    s = [f'<svg class="fig" viewBox="0 0 {W} {H}" role="img" aria-label="'
         f'Relative excess after an SAA exit, by energy channel">']
    s.append(f'<rect class="frame" x="{M["l"]}" y="{M["t"]}" width="{iw}" '
             f'height="{ih}" />')
    s.append('<g class="axis">')
    for v in axis_ticks(0, max(t), 6):
        if v > max(t):
            continue
        s.append(f'<line x1="{px(v):.1f}" y1="{M["t"]}" x2="{px(v):.1f}" '
                 f'y2="{M["t"]+ih}"/>')
        s.append(f'<text x="{px(v):.1f}" y="{M["t"]+ih+14}" text-anchor="middle">'
                 f'{v:g}</text>')
    for v in axis_ticks(y0, y1, 4):
        if not (y0 <= v <= y1):
            continue
        s.append(f'<line x1="{M["l"]}" y1="{py(v):.1f}" x2="{M["l"]+iw}" '
                 f'y2="{py(v):.1f}"/>')
        s.append(f'<text x="{M["l"]-7}" y="{py(v)+3:.1f}" text-anchor="end">'
                 f'{v:g}</text>')
    s.append(f'<text x="{M["l"]+iw/2}" y="{H-4}" text-anchor="middle">'
             f'minutes since leaving the SAA</text>')
    s.append(f'<text x="13" y="{M["t"]+ih/2}" text-anchor="middle" '
             f'transform="rotate(-90 13 {M["t"]+ih/2})">rate / quiet rate</text>')
    s.append('</g>')
    if y0 < 1.0 < y1:
        s.append(f'<line class="unity" x1="{M["l"]}" y1="{py(1.0):.1f}" '
                 f'x2="{M["l"]+iw}" y2="{py(1.0):.1f}" />')
    wanted = []
    for i, row in enumerate(rows):
        frac = i / max(len(rows) - 1, 1)
        # Soft channels blue, hard channels red, matching the viewer.
        col = (f"rgb({int(47+170*frac)},{int(111-28*frac)},{int(208-129*frac)})")
        s.append(f'<path class="line" style="stroke:{col}" '
                 f'd="{path_from(t, row["rel"], px, py)}" />')
        last = next((v for v in reversed(row["rel"]) if v is not None), None)
        if last is not None:
            wanted.append((py(last),
                           (f'{row["lo"]:g}-{row["hi"]:g} keV', col)))
    for y, (text, col) in place_labels(wanted, M["t"] + 11, M["t"] + ih - 2):
        s.append(f'<text class="serieslab" x="{M["l"]+iw+6}" y="{y+3:.1f}" '
                 f'style="fill:{col}">{text}</text>')
    s.append('</svg>')
    return "".join(s)


def main():
    with open(SRC) as handle:
        data = json.load(handle)

    fit = data["fit"]
    nai_tau = sorted(fit[d]["tau"] for d in NAI if d in fit)
    rows = []
    for det in NAI + BGO:
        f = fit.get(det)
        if not f:
            continue
        rows.append(
            f'<tr{" class=bgo" if det in BGO else ""}><td class="m">{det}</td>'
            f'<td class="n">{f["A"]:,.0f}</td><td class="n">{f["tau"]:.2f}</td>'
            f'<td class="n">{f["B"]:,.0f}</td>'
            f'<td class="n">{f["A"]/f["B"]:.2f}</td>'
            f'<td class="n">{f["peak_excess"]:+.1f}%</td></tr>')

    # Excess per energy channel at a few delays, for the table under the chart.
    import numpy as np
    t_arr = np.asarray(data["t"], dtype=float)
    chan_rows = []
    for row in data["chan"]["rows"]:
        rel = np.array([np.nan if v is None else v for v in row["rel"]])

        def excess(lo, hi):
            m = (t_arr >= lo) & (t_arr < hi)
            v = np.nanmedian(rel[m]) if m.any() else np.nan
            return "&mdash;" if not np.isfinite(v) else f"{100*(v-1):+.1f}%"

        chan_rows.append(
            f'<tr><td class="m">{row["lo"]:g}&ndash;{row["hi"]:g} keV</td>'
            f'<td class="n">{row["base"]:,.0f}</td>'
            f'<td class="n">{excess(0,1)}</td>'
            f'<td class="n">{excess(2,5)}</td>'
            f'<td class="n">{excess(10,20)}</td></tr>')
    chan_rows = "".join(chan_rows)

    day = data["day"]
    breaks = sum(1 for v in day["r"] if v is None)
    n_saa = len(day["saa"])
    covered = sum(1 for v in day["r"] if v is not None) * day["step"] / 3600.0

    built = datetime.date.today().isoformat()
    html = f"""<title>SAA Background in Fermi GBM</title>
<style>
  :root {{
    --paper:#f5f7fa; --surface:#ffffff; --surface-2:#eef1f6;
    --ink:#141b26; --ink-soft:#46525f; --muted:#6b7686; --rule:#d8dee8;
    --accent:{ACCENT}; --hot:{HOT}; --warm:{WARM}; --cool:{COOL};
    --grid:#c9d2de;
    --shadow:0 1px 2px rgba(20,27,38,.06), 0 6px 18px rgba(20,27,38,.05);
    --serif: ui-serif, Georgia, "Iowan Old Style", "Times New Roman", serif;
    --sans: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
    --mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
      --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
      --accent:#6fa3ea; --hot:#ef8079; --warm:#ffd166; --cool:#4cc79a;
      --grid:#2b3747;
      --shadow:0 1px 2px rgba(0,0,0,.4), 0 8px 22px rgba(0,0,0,.35);
    }}
  }}
  :root[data-theme="dark"] {{
    --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
    --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
    --accent:#6fa3ea; --hot:#ef8079; --warm:#ffd166; --cool:#4cc79a;
    --grid:#2b3747;
  }}
  body {{ margin:0; background:var(--paper); color:var(--ink);
    font-family:var(--serif); line-height:1.6; }}
  .wrap {{ max-width:1100px; margin:0 auto; padding:2.2rem 1.3rem 4rem; }}
  h1 {{ font-size:1.85rem; line-height:1.2; margin:0 0 .3rem; }}
  h2 {{ font-size:1.2rem; margin:2.2rem 0 .5rem; font-family:var(--sans); }}
  .sub {{ color:var(--muted); font-size:.92rem; margin:0 0 1.6rem;
    font-family:var(--sans); }}
  p {{ max-width:74ch; }}
  .card {{ background:var(--surface); border:1px solid var(--rule);
    border-radius:12px; padding:1rem 1.1rem; box-shadow:var(--shadow);
    margin:1rem 0 1.4rem; }}
  svg.fig {{ width:100%; height:auto; display:block; }}
  svg.fig text {{ font-family:var(--sans); font-size:10px; fill:var(--muted); }}
  svg.fig .frame {{ fill:none; stroke:var(--rule); stroke-width:1; }}
  svg.fig .axis line {{ stroke:var(--grid); stroke-width:.5; }}
  svg.fig .line {{ fill:none; stroke:var(--accent); stroke-width:1.4;
    stroke-linejoin:round; }}
  svg.fig .saaband {{ fill:var(--hot); fill-opacity:.17; }}
  svg.fig .unity {{ stroke:var(--muted); stroke-width:1;
    stroke-dasharray:4,3; opacity:.7; }}
  svg.fig .serieslab {{ font-family:var(--mono); font-size:9.5px; }}
  .cap {{ font-size:.78rem; color:var(--muted); font-family:var(--sans);
    margin:.5rem 0 0; max-width:86ch; }}
  table {{ border-collapse:collapse; font-family:var(--sans);
    font-size:.86rem; width:100%; }}
  th, td {{ padding:.34rem .6rem; border-bottom:1px solid var(--rule);
    text-align:left; }}
  th {{ font-size:.74rem; text-transform:uppercase; letter-spacing:.04em;
    color:var(--muted); font-weight:650; }}
  td.n {{ text-align:right; font-family:var(--mono); }}
  td.m {{ font-family:var(--mono); font-weight:650; }}
  tr.bgo td {{ color:var(--hot); }}
  .tiles {{ display:grid; grid-template-columns:repeat(auto-fit,minmax(190px,1fr));
    gap:.8rem; margin:1.2rem 0; }}
  .tile {{ background:var(--surface-2); border-radius:10px; padding:.8rem .9rem; }}
  .tile .v {{ display:block; font-size:1.5rem; font-weight:680;
    font-family:var(--sans); }}
  .tile .k {{ display:block; font-size:.76rem; color:var(--muted);
    font-family:var(--sans); margin-top:.15rem; }}
  .warn {{ border-left:3px solid var(--warm); padding-left:.9rem;
    background:var(--surface-2); border-radius:0 8px 8px 0; }}
  footer {{ margin-top:3rem; color:var(--muted); font-size:.8rem;
    font-family:var(--sans); border-top:1px solid var(--rule);
    padding-top:1rem; }}
  a {{ color:var(--accent); }}
</style>

<div class="wrap">
<h1>Background from the South Atlantic Anomaly</h1>
<p class="sub">{data['days']} days of continuous GBM data, {data['first']} to
  {data['last']} &middot; {data['exits']} SAA passages &middot; built {built}</p>

<p>Fermi passes through the South Atlantic Anomaly several times a day. The
  detectors are switched off while it does, and when they come back on they
  are counting their own housing: charged particles have activated the
  spacecraft, and that activation decays away over the following minutes. It
  is a real background component, and it is one of the things a physical
  background model has to get right.</p>

<p>None of this is visible anywhere else in this project. Every other
  measurement here is trigger-centred, covering a few hundred seconds around
  a burst, and the activation decays over a span many times longer than that.
  Inside such a window it is an almost constant offset, absorbed whole by the
  background polynomial's constant term. Seeing it at all needs continuous
  days.</p>

<div class="tiles">
  <div class="tile"><span class="v">{data['exits']}</span>
    <span class="k">SAA passages measured</span></div>
  <div class="tile"><span class="v">+{fit['n0']['peak_excess']:.0f}%</span>
    <span class="k">excess on n0 right after an exit</span></div>
  <div class="tile"><span class="v">{nai_tau[len(nai_tau)//2]:.1f} min</span>
    <span class="k">median NaI decay constant</span></div>
  <div class="tile"><span class="v">{fit['b0']['tau']:.0f} min</span>
    <span class="k">BGO decay constant, b0</span></div>
</div>

<h2>One day, as it happens</h2>
<div class="card">
  {chart_day(data['day'])}
  <p class="cap">Detector {data['day']['det']} across {data['day']['day']},
    binned to {data['day']['step']:.0f} s. Shaded bands are SAA passages, and
    the trace genuinely stops inside each one because the detectors are
    switched off &mdash; {breaks} breaks for {n_saa} passages, leaving
    {covered:.1f} of the 24 hours with data. Each restart is followed by a
    spike that decays within minutes. The slower rise and fall underneath is
    orbital.</p>
</div>

<h2>Stacked on the moment of exit</h2>
<div class="card">
  {chart_fold(data, NAI, 'var(--accent)')}
  <p class="cap">All twelve NaI detectors, every passage stacked on the
    instant Fermi left the SAA. Median rate, not mean: the quiet reference
    window carries the most exposure and therefore the most bursts and
    flares, and taking means put that supposed floor <em>above</em> the
    profile it was meant to anchor.</p>
</div>

<div class="card">
  {chart_fold(data, BGO, 'var(--hot)')}
  <p class="cap">The two BGOs on the same axes. They behave differently
    enough that sharing one decay constant across all fourteen detectors
    would be wrong.</p>
</div>

<h2>Is the activation harder than the background?</h2>
<p>Markedly so. Splitting the same fold by CTIME's eight energy channels, the
  excess in the first minute climbs from about a tenth in the softest channel
  to roughly three quarters around half an MeV &mdash; where the de-excitation
  lines of activated material sit. The two softest channels go <em>below</em>
  their quiet level after a few minutes, which is orbital rather than
  activation: an SAA exit happens at a particular point in the orbit, and the
  soft background tracks orbital position.</p>
<div class="card">
  {chart_chan(data['chan'], data['t'])}
  <p class="cap">Detector {data['chan']['det']}, each channel divided by its
    own quiet level so channels differing by orders of magnitude in absolute
    rate can be compared. Dashed line is the quiet level. Soft channels blue,
    hard channels red.</p>
</div>
<div class="card">
<table>
  <tr><th>Channel</th><th style="text-align:right">Quiet rate</th>
      <th style="text-align:right">0&ndash;1 min</th>
      <th style="text-align:right">2&ndash;5 min</th>
      <th style="text-align:right">10&ndash;20 min</th></tr>
  {chan_rows}
</table>
  <p class="cap">Excess over each channel's own quiet level.</p>
</div>

<h2>Fitted decay</h2>
<div class="card">
<table>
  <tr><th>Detector</th><th style="text-align:right">Amplitude A</th>
      <th style="text-align:right">tau (min)</th>
      <th style="text-align:right">Baseline B</th>
      <th style="text-align:right">A/B</th>
      <th style="text-align:right">Peak excess</th></tr>
  {''.join(rows)}
</table>
  <p class="cap">A single exponential A&middot;exp(&minus;t/tau) + B fitted
    over the first hour. Peak excess is the first half-minute measured against
    the 60&ndash;90 minute level.</p>
</div>

<h2>What this is not</h2>
<div class="card warn">
  <p>A single exponential is a convenience rather than the truth. The profile
    dips about 4% <em>below</em> the 60&ndash;90 minute reference around
    30&ndash;45 minutes, so the baseline is not flat &mdash; almost certainly
    orbital, since geomagnetic latitude varies on the same roughly 96 minute
    period and is itself a background driver. The fit also only sees the fast
    component; a slower one takes hours and needs a whole day modelled at
    once.</p>
  <p>These numbers are a direct measurement, not a model. Nothing here is
    fitted jointly with the other background components &mdash; cosmic
    diffuse, Earth albedo, point sources &mdash; which is what a physical
    background model such as
    <a href="https://arxiv.org/abs/2005.11219">Biltzinger et al. (2020)</a>
    does, and which is the obvious next step.</p>
</div>

<footer>
  Fermi GBM daily CTIME and position history from the
  <a href="https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/">HEASARC
  archive</a>. Part of a wider comparison of reproduced and official GBM
  burst localizations.
</footer>
</div>
"""
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as handle:
        handle.write(html)
    non_ascii = sum(1 for ch in html if ord(ch) > 127)
    print(f"wrote {OUT} ({len(html)/1024:.0f} KB, non-ascii={non_ascii})")


if __name__ == "__main__":
    main()
