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
    """Pick a round step giving close to n intervals.

    The obvious rule -- take the first nice step at or above (hi-lo)/n --
    overshoots badly. For a span of 96 asking for 3 intervals it rejects 25
    as too small and jumps to 50, so the panel ends up with two gridlines.
    Scoring every candidate on how near it lands to the requested count, and
    breaking ties toward more ticks, keeps the axis readable.
    """
    import math
    span = hi - lo
    if span <= 0:
        return [lo]
    def ticks_for(step):
        out, v = [], math.ceil(lo / step) * step
        while v <= hi + 1e-9:
            out.append(round(v, 6))
            v += step
        return out

    mag = 10 ** math.floor(math.log10(span / max(n, 1)))
    best, best_score = None, None
    for mult in (0.5, 1, 2, 2.5, 5, 10, 20):
        step = mult * mag
        got = ticks_for(step)
        # Score the ticks that actually land inside the range, not the
        # nominal interval count: a step dividing the span 2.6 times can
        # still place only two labels once it is aligned to round numbers,
        # which is how b0 ended up with a single gridline.
        if not (3 <= len(got) <= 9):
            continue
        score = (abs(len(got) - (n + 1)), -len(got))
        if best_score is None or score < best_score:
            best, best_score = got, score
    return best if best is not None else ticks_for(span / max(n, 1))


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


def chan_colour(frac):
    """Soft blue to hard red, matching the viewer's band colours."""
    return f"rgb({int(47+170*frac)},{int(111-28*frac)},{int(208-129*frac)})"


def grid_panel(det, rows, t, y0, y1, W=236, H=150):
    """One small multiple: every energy channel of one detector."""
    M = dict(l=34, r=6, t=16, b=22)
    iw, ih = W - M["l"] - M["r"], H - M["t"] - M["b"]
    px = lambda v: M["l"] + iw * v / (max(t) or 1)
    py = lambda v: M["t"] + ih * (y1 - v) / (y1 - y0 or 1)
    s = [f'<svg class="mini" viewBox="0 0 {W} {H}" role="img" '
         f'aria-label="Detector {det}, excess after an SAA exit by energy">']
    s.append(f'<rect class="frame" x="{M["l"]}" y="{M["t"]}" width="{iw}" '
             f'height="{ih}" />')
    s.append('<g class="axis">')
    for v in (0, 30, 60, 90):
        if v > max(t):
            continue
        s.append(f'<line x1="{px(v):.1f}" y1="{M["t"]}" x2="{px(v):.1f}" '
                 f'y2="{M["t"]+ih}"/>')
        s.append(f'<text x="{px(v):.1f}" y="{M["t"]+ih+12}" '
                 f'text-anchor="middle">{v}</text>')
    for v in axis_ticks(y0, y1, 3):
        if not (y0 <= v <= y1):
            continue
        s.append(f'<line x1="{M["l"]}" y1="{py(v):.1f}" x2="{M["l"]+iw}" '
                 f'y2="{py(v):.1f}"/>')
        s.append(f'<text x="{M["l"]-5}" y="{py(v)+3:.1f}" text-anchor="end">'
                 f'{v:g}</text>')
    s.append('</g>')
    if y0 < 1.0 < y1:
        s.append(f'<line class="unity" x1="{M["l"]}" y1="{py(1.0):.1f}" '
                 f'x2="{M["l"]+iw}" y2="{py(1.0):.1f}" />')
    for i, row in enumerate(rows):
        frac = i / max(len(rows) - 1, 1)
        s.append(f'<path class="line thin" style="stroke:{chan_colour(frac)}" '
                 f'd="{path_from(t, row["rel"], px, py)}" />')
    s.append(f'<text class="minilab" x="{M["l"]}" y="11">{det}</text>')
    s.append(f'<text class="minisub" x="{W-M["r"]}" y="11" text-anchor="end">'
             f'{rows[0]["lo"]:g}&#8211;{rows[-1]["hi"]:g} keV</text>')
    s.append('</svg>')
    return "".join(s)


def chan_legend(rows, width=980):
    """Shared key for the grid: which colour is which energy."""
    s = [f'<svg class="key" viewBox="0 0 {width} 26" role="img" '
         f'aria-label="Energy channel colour key">']
    x = 4
    for i, row in enumerate(rows):
        frac = i / max(len(rows) - 1, 1)
        s.append(f'<line x1="{x}" y1="13" x2="{x+20}" y2="13" '
                 f'style="stroke:{chan_colour(frac)};stroke-width:2.4" />')
        label = f'{row["lo"]:g}-{row["hi"]:g}'
        s.append(f'<text x="{x+25}" y="16">{label}</text>')
        x += 32 + 6.1 * len(label)
    s.append(f'<text x="{min(x+6, width-4)}" y="16">keV</text>')
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

    # The all-detector grid. One shared vertical scale across the twelve NaI
    # so panels can be compared at a glance; the BGOs get their own, because
    # their excess is roughly twice as large and would flatten the rest.
    chan_all = data.get("chan_all") or {}

    def span(dets):
        vals = [v for d in dets for row in chan_all.get(d, [])
                for v in row["rel"] if v is not None]
        if not vals:
            return 0.9, 1.1
        return min(vals) - 0.02, max(vals) + 0.02

    # The explorer draws client-side, so it needs the folds rather than
    # pre-rendered SVG. Only the fields it actually reads are shipped.
    explore_json = json.dumps({
        "t": data["t"],
        "nai": [d for d in NAI if d in chan_all],
        "bgo": [d for d in BGO if d in chan_all],
        "chan": {d: [{"lo": r["lo"], "hi": r["hi"], "base": r["base"],
                      "rel": r["rel"]} for r in rows]
                 for d, rows in chan_all.items()},
    }, separators=(",", ":"))

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
  .grid {{ display:grid; gap:.4rem;
    grid-template-columns:repeat(4,1fr); }}
  @media (max-width:860px) {{ .grid {{ grid-template-columns:repeat(3,1fr); }} }}
  @media (max-width:620px) {{ .grid {{ grid-template-columns:repeat(2,1fr); }} }}
  svg.mini {{ width:100%; height:auto; display:block; }}
  svg.mini text {{ font-family:var(--sans); font-size:8.5px; fill:var(--muted); }}
  svg.mini .frame {{ fill:none; stroke:var(--rule); stroke-width:1; }}
  svg.mini .axis line {{ stroke:var(--grid); stroke-width:.5; }}
  svg.mini .line.thin {{ fill:none; stroke-width:1; stroke-linejoin:round; }}
  svg.mini .unity {{ stroke:var(--muted); stroke-width:.8;
    stroke-dasharray:3,2; opacity:.7; }}
  svg.mini .minilab {{ font-family:var(--mono); font-size:10px;
    font-weight:650; fill:var(--ink); }}
  svg.mini .minisub {{ font-size:8px; opacity:.8; }}
  svg.key {{ width:100%; height:auto; display:block; margin:0 0 .5rem; }}
  svg.key text {{ font-family:var(--mono); font-size:9px; fill:var(--muted); }}

  /* ------------------------------------------------------- controls */
  .ctl {{ display:flex; flex-direction:column; gap:.5rem;
    padding-bottom:.7rem; margin-bottom:.4rem;
    border-bottom:1px solid var(--rule); }}
  .ctlrow {{ display:flex; flex-wrap:wrap; align-items:center; gap:.3rem; }}
  .ctllab {{ font-family:var(--sans); font-size:.72rem; font-weight:650;
    text-transform:uppercase; letter-spacing:.04em; color:var(--muted);
    min-width:68px; }}
  .chips {{ display:flex; flex-wrap:wrap; gap:.25rem; }}
  .chips button, .quick button {{ font-family:var(--mono); font-size:.76rem;
    padding:.2rem .5rem; border-radius:6px; cursor:pointer;
    border:1px solid var(--rule); background:var(--surface-2);
    color:var(--ink-soft); }}
  .chips button[aria-pressed="true"] {{ background:var(--accent);
    border-color:var(--accent); color:#fff; font-weight:650; }}
  .quick {{ display:flex; gap:.25rem; margin-left:.5rem; }}
  .quick button {{ font-family:var(--sans); font-size:.72rem;
    background:transparent; color:var(--muted); }}
  .quick button:hover, .chips button:hover {{ border-color:var(--accent); }}
  /* Energy chips carry their channel's own colour when selected. */
  .chips button.ch[aria-pressed="true"] {{ color:#fff; font-weight:650; }}
  #explore:empty::after {{ content:"Nothing selected."; color:var(--muted);
    font-family:var(--sans); font-size:.82rem; }}
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

<h2>Every detector, every energy range</h2>
<p>Pick the detectors and the energy channels. Each detector gets its own
  panel. Nothing is drawn for the time Fermi spends inside the SAA: the
  detectors are off, so there is no measurement to show.</p>

<div class="card">
  <div class="ctl">
    <div class="ctlrow">
      <span class="ctllab">Detectors</span>
      <span id="detchips" class="chips"></span>
      <span class="quick"><button data-pick="nai">all NaI</button><button
        data-pick="bgo">both BGO</button><button data-pick="few">just n0 n6 b0</button><button
        data-pick="none">clear</button></span>
    </div>
    <div class="ctlrow">
      <span class="ctllab">Energy</span>
      <span id="chanchips" class="chips"></span>
      <span class="quick"><button data-ch="all">all</button><button
        data-ch="hard">hard only</button><button data-ch="soft">soft only</button></span>
    </div>
    <div class="ctlrow">
      <span class="ctllab">Window</span>
      <span class="chips" id="spanchips"></span>
      <span class="ctllab" style="margin-left:1.2rem">Vertical</span>
      <span class="chips" id="modechips"></span>
      <span class="ctllab" style="margin-left:1.2rem">Scale</span>
      <span class="chips" id="scalechips"></span>
    </div>
  </div>
  <p id="explorenote" class="cap"></p>
  <div id="explore" class="grid"></div>
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

<script>
(function () {{
  var D = {explore_json};
  var NAI = D.nai, BGO = D.bgo, ALL = NAI.concat(BGO);
  var state = {{
    dets: ['n0','n3','n6','n9','b0'],
    chans: [0,1,2,3,4,5,6,7],
    span: 30, mode: 'pct', scale: 'shared'
  }};

  function colour(i) {{
    var f = i / 7;
    return 'rgb(' + Math.round(47+170*f) + ',' + Math.round(111-28*f) +
           ',' + Math.round(208-129*f) + ')';
  }}
  // Channel edges differ between NaI and BGO, so a chip shows the NaI range
  // and each panel states its own underneath.
  function chanLabel(i) {{
    var r = D.chan['n0'][i];
    return (r.lo < 1000 ? r.lo.toFixed(0) : (r.lo/1000).toFixed(0) + 'k') + '-' +
           (r.hi < 1000 ? r.hi.toFixed(0) : (r.hi/1000).toFixed(0) + 'k');
  }}

  function chip(host, label, pressed, onclick, cls, style) {{
    var b = document.createElement('button');
    b.textContent = label;
    b.setAttribute('aria-pressed', pressed ? 'true' : 'false');
    if (cls) b.className = cls;
    if (style) b.setAttribute('style', style);
    b.addEventListener('click', onclick);
    host.appendChild(b);
    return b;
  }}

  function buildControls() {{
    var dh = document.getElementById('detchips');
    dh.innerHTML = '';
    ALL.forEach(function (d) {{
      chip(dh, d, state.dets.indexOf(d) >= 0, function () {{
        var i = state.dets.indexOf(d);
        if (i >= 0) state.dets.splice(i, 1); else state.dets.push(d);
        state.dets.sort(function (a, b) {{ return ALL.indexOf(a) - ALL.indexOf(b); }});
        render();
      }});
    }});
    var ch = document.getElementById('chanchips');
    ch.innerHTML = '';
    for (var i = 0; i < 8; i++) {{
      (function (i) {{
        var on = state.chans.indexOf(i) >= 0;
        chip(ch, chanLabel(i), on, function () {{
          var k = state.chans.indexOf(i);
          if (k >= 0) state.chans.splice(k, 1); else state.chans.push(i);
          state.chans.sort(function (a, b) {{ return a - b; }});
          render();
        }}, 'ch', on ? 'background:' + colour(i) + ';border-color:' + colour(i) : '');
      }})(i);
    }}
    var sh = document.getElementById('spanchips');
    sh.innerHTML = '';
    [10, 30, 90].forEach(function (v) {{
      chip(sh, v + ' min', state.span === v, function () {{
        state.span = v; render();
      }});
    }});
    var mh = document.getElementById('modechips');
    mh.innerHTML = '';
    [['pct', '% above quiet'], ['abs', 'counts/s']].forEach(function (m) {{
      chip(mh, m[1], state.mode === m[0], function () {{
        state.mode = m[0]; render();
      }});
    }});
    var kh = document.getElementById('scalechips');
    kh.innerHTML = '';
    [['shared', 'shared'], ['auto', 'per panel']].forEach(function (m) {{
      chip(kh, m[1], state.scale === m[0], function () {{
        state.scale = m[0]; render();
      }});
    }});
  }}

  document.querySelectorAll('.quick button[data-pick]').forEach(function (b) {{
    b.addEventListener('click', function () {{
      var p = b.getAttribute('data-pick');
      state.dets = p === 'nai' ? NAI.slice() : p === 'bgo' ? BGO.slice()
                 : p === 'few' ? ['n0','n6','b0'] : [];
      render();
    }});
  }});
  document.querySelectorAll('.quick button[data-ch]').forEach(function (b) {{
    b.addEventListener('click', function () {{
      var p = b.getAttribute('data-ch');
      state.chans = p === 'all' ? [0,1,2,3,4,5,6,7]
                  : p === 'hard' ? [4,5,6,7] : [0,1,2,3];
      render();
    }});
  }});

  function valueAt(det, c, i) {{
    var row = D.chan[det][c], v = row.rel[i];
    if (v === null) return null;
    return state.mode === 'pct' ? (v - 1) * 100 : v * row.base;
  }}

  function panel(det, y0, y1, w, h) {{
    var M = {{l: 40, r: 7, t: 17, b: 23}};
    var iw = w - M.l - M.r, ih = h - M.t - M.b;
    var t = D.t, span = state.span;
    var px = function (v) {{ return M.l + iw * v / span; }};
    var py = function (v) {{ return M.t + ih * (y1 - v) / (y1 - y0 || 1); }};
    var s = ['<svg class="mini" viewBox="0 0 ' + w + ' ' + h + '" role="img" ' +
             'aria-label="Detector ' + det + ' after an SAA exit">'];
    s.push('<rect class="frame" x="' + M.l + '" y="' + M.t + '" width="' + iw +
           '" height="' + ih + '" />');
    s.push('<g class="axis">');
    var step = span <= 10 ? 2 : span <= 30 ? 10 : 30;
    for (var v = 0; v <= span; v += step) {{
      s.push('<line x1="' + px(v).toFixed(1) + '" y1="' + M.t + '" x2="' +
             px(v).toFixed(1) + '" y2="' + (M.t + ih) + '"/>');
      s.push('<text x="' + px(v).toFixed(1) + '" y="' + (M.t + ih + 12) +
             '" text-anchor="middle">' + v + '</text>');
    }}
    var ticks = niceTicks(y0, y1, 3);
    ticks.forEach(function (v) {{
      s.push('<line x1="' + M.l + '" y1="' + py(v).toFixed(1) + '" x2="' +
             (M.l + iw) + '" y2="' + py(v).toFixed(1) + '"/>');
      s.push('<text x="' + (M.l - 5) + '" y="' + (py(v) + 3).toFixed(1) +
             '" text-anchor="end">' + fmt(v) + '</text>');
    }});
    s.push('</g>');
    var zero = state.mode === 'pct' ? 0 : null;
    if (zero !== null && y0 < zero && zero < y1) {{
      s.push('<line class="unity" x1="' + M.l + '" y1="' + py(zero).toFixed(1) +
             '" x2="' + (M.l + iw) + '" y2="' + py(zero).toFixed(1) + '" />');
    }}
    state.chans.forEach(function (c) {{
      var d = '', pen = false;
      for (var i = 0; i < t.length; i++) {{
        if (t[i] > span) break;
        var v = valueAt(det, c, i);
        if (v === null) {{ pen = false; continue; }}
        d += (pen ? 'L' : 'M') + px(t[i]).toFixed(1) + ',' + py(v).toFixed(1);
        pen = true;
      }}
      if (d) s.push('<path class="line thin" style="stroke:' + colour(c) +
                    '" d="' + d + '" />');
    }});
    s.push('<text class="minilab" x="' + M.l + '" y="12">' + det + '</text>');
    var rows = D.chan[det], sel = state.chans;
    if (sel.length) {{
      var lo = rows[Math.min.apply(null, sel)].lo;
      var hi = rows[Math.max.apply(null, sel)].hi;
      s.push('<text class="minisub" x="' + (w - M.r) + '" y="12" ' +
             'text-anchor="end">' + fmtkev(lo) + '-' + fmtkev(hi) + ' keV</text>');
    }}
    s.push('</svg>');
    return s.join('');
  }}

  function fmtkev(v) {{
    return v >= 1000 ? (v / 1000).toFixed(v >= 10000 ? 0 : 1) + 'k' : v.toFixed(0);
  }}
  function fmt(v) {{
    var r = Math.round(v * 10) / 10;
    return Number.isInteger(r) ? String(r) : r.toFixed(1);
  }}
  // Same scoring as the Python axis_ticks: taking the first nice step at or
  // above the raw interval overshoots and leaves two gridlines on a panel.
  function niceTicks(lo, hi, n) {{
    var span = hi - lo;
    if (!(span > 0)) return [lo];
    function ticksFor(step) {{
      var out = [], v = Math.ceil(lo / step) * step;
      for (; v <= hi + 1e-9; v += step) out.push(Math.round(v * 1e6) / 1e6);
      return out;
    }}
    var mag = Math.pow(10, Math.floor(Math.log(span / n) / Math.LN10));
    var best = null, bestScore = null;
    [0.5, 1, 2, 2.5, 5, 10, 20].forEach(function (m) {{
      var got = ticksFor(m * mag);
      if (got.length < 3 || got.length > 9) return;
      var score = Math.abs(got.length - (n + 1)) - got.length / 1000;
      if (bestScore === null || score < bestScore) {{ best = got; bestScore = score; }}
    }});
    return best === null ? ticksFor(span / n) : best;
  }}

  function range(dets) {{
    var lo = Infinity, hi = -Infinity;
    dets.forEach(function (det) {{
      state.chans.forEach(function (c) {{
        for (var i = 0; i < D.t.length; i++) {{
          if (D.t[i] > state.span) break;
          var v = valueAt(det, c, i);
          if (v === null) continue;
          if (v < lo) lo = v;
          if (v > hi) hi = v;
        }}
      }});
    }});
    if (!isFinite(lo)) return [0, 1];
    var pad = (hi - lo) * 0.08 || Math.abs(hi) * 0.08 || 1;
    return [lo - pad, hi + pad];
  }}

  function render() {{
    buildControls();
    var host = document.getElementById('explore');
    host.innerHTML = '';
    var note = document.getElementById('explorenote');
    if (!state.dets.length || !state.chans.length) {{
      note.textContent = 'Pick at least one detector and one energy channel.';
      return;
    }}
    // Mixing NaI and BGO on one scale is misleading -- their channel edges
    // and their excess are both quite different -- so a shared scale is
    // shared within a type, not across.
    var nai = state.dets.filter(function (d) {{ return NAI.indexOf(d) >= 0; }});
    var bgo = state.dets.filter(function (d) {{ return BGO.indexOf(d) >= 0; }});
    var rn = range(nai), rb = range(bgo);
    var n = state.dets.length;
    var cols = n <= 1 ? 1 : n <= 2 ? 2 : n <= 6 ? 3 : 4;
    host.style.gridTemplateColumns = 'repeat(' + cols + ',1fr)';
    var w = 236, h = cols <= 2 ? 190 : 150;
    state.dets.forEach(function (det) {{
      var isB = BGO.indexOf(det) >= 0;
      var r = state.scale === 'auto' ? range([det]) : (isB ? rb : rn);
      host.insertAdjacentHTML('beforeend', panel(det, r[0], r[1], w, h));
    }});
    var unit = state.mode === 'pct' ? 'per cent above the quiet level'
                                    : 'counts/s';
    note.innerHTML = state.dets.length + ' detector' +
      (state.dets.length > 1 ? 's' : '') + ', ' + state.chans.length +
      ' channel' + (state.chans.length > 1 ? 's' : '') + ', first ' +
      state.span + ' minutes after the exit, vertical axis in ' + unit +
      ', ' + (state.scale === 'auto' ? 'each panel scaled to its own range'
              : 'one scale shared across panels of the same detector type') +
      '. Soft channels blue, hard red.';
  }}

  render();
}})();
</script>

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
