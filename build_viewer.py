"""
Build the single-page burst viewer: pick a class, pick a GRB, see its plot.

Every one of the ~2,000 bursts is selectable. That is only possible because
the sky regions are stored as contour polygons and drawn as SVG in the
browser; the PNG gallery it replaces totalled 81 MB against a 16 MB page
budget, so it could only ever show a few hundred.

Output: site/viewer.html
"""
import json
import glob
import os
import statistics
import datetime

SRC = "site/vec/bursts"
OUT = "site/viewer.html"

OFFICIAL_C = "#2f6fd0"
T90_C = "#e8913a"
PEAK_C = "#d9534f"
EXT_C = "#1f9e6e"


def main():
    records = []
    for path in sorted(glob.glob(os.path.join(SRC, "*.json"))):
        with open(path) as handle:
            records.append(json.load(handle))
    if not records:
        raise SystemExit("no contour data; run extract_contours.py first")

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
    --good:#2e7d55; --crit:#bf3b36; --grid:#c9d2de;
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
      --good:#58b98a; --crit:#e8706a; --grid:#2f3b4d;
      --shadow:0 1px 2px rgba(0,0,0,.4), 0 6px 18px rgba(0,0,0,.3);
    }}
  }}
  :root[data-theme="dark"] {{
    --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
    --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
    --accent:#6fa3ea; --t90:#f0a85f; --peak:#ef8079; --ext:#4cc79a;
    --good:#58b98a; --crit:#e8706a; --grid:#2f3b4d;
    --shadow:0 1px 2px rgba(0,0,0,.4), 0 6px 18px rgba(0,0,0,.3);
  }}
  *,*::before,*::after {{ box-sizing:border-box; }}
  body {{ margin:0; background:var(--paper); color:var(--ink);
    font-family:var(--sans); font-size:16px; line-height:1.6;
    -webkit-font-smoothing:antialiased; }}
  .wrap {{ max-width:1080px; margin:0 auto;
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
  .stage {{ display:grid; grid-template-columns:minmax(0,1.15fr) minmax(0,1fr);
    gap:1.4rem; align-items:start; }}
  @media (max-width:820px) {{ .stage {{ grid-template-columns:1fr; }} }}
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
  .star {{ color:var(--ext); font-size:1rem; line-height:1; }}

  /* ---------------------------------------------------------- facts */
  .facts {{ background:var(--surface); border:1px solid var(--rule);
    border-radius:12px; padding:1.1rem 1.2rem; box-shadow:var(--shadow); }}
  .facts h3 {{ margin:0 0 .2rem; font-size:1.05rem; font-family:var(--mono); }}
  .facts .sub {{ color:var(--muted); font-size:.82rem; margin:0 0 .9rem; }}
  dl.kv {{ margin:0; display:grid; grid-template-columns:auto 1fr; gap:.3rem .9rem; }}
  dl.kv dt {{ color:var(--muted); font-size:.82rem; }}
  dl.kv dd {{ margin:0; font-size:.84rem; text-align:right;
    font-variant-numeric:tabular-nums; }}
  .sep {{ border-top:1px solid var(--rule); margin:.9rem 0; }}
  .big {{ font-family:var(--serif); font-size:1.5rem; line-height:1.1;
    font-variant-numeric:tabular-nums; }}
  .big.good {{ color:var(--good); }} .big.bad {{ color:var(--crit); }}
  .delta {{ font-size:.82rem; }}
  .delta.good {{ color:var(--good); }} .delta.bad {{ color:var(--crit); }}
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
      <svg class="sky" id="sky" viewBox="0 0 440 400" role="img"
           aria-label="Sky localization comparison for the selected burst"></svg>
      <div class="legend">
        <span class="key"><span class="sw off"></span> official Fermi GBM 90%</span>
        <span class="key"><span class="sw t90"></span> mine, full T90</span>
        <span class="key"><span class="sw pk"></span> mine, best S/N</span>
        <span class="key"><span class="star">&#9733;</span> other mission</span>
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
      ext: v('--ext', '{EXT_C}')
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
      var tag = r.x ? ' *' : '';
      return '<option value="' + r.t + '">' + r.t + '  (' + r.ps.toFixed(2) +
             '\\u00B0' + tag + ')</option>';
    }}).join('');
    if (keep && list.some(function (r) {{ return r.t === keep; }})) grbSel.value = keep;
    pickcount.textContent = list.length + ' bursts  \\u00B7  * = has external position';
  }}

  // ---------------------------------------------------------------- plot
  function wrapNear(v, c) {{ return c + ((v - c + 180) % 360 + 360) % 360 - 180; }}

  function draw(r) {{
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
    h += '<h3>' + esc(r.t) + '</h3>';
    h += '<p class="sub">' + esc(r.n) + ' &middot; T90 ' + r.t90 + ' s (' + r.k +
         ') &middot; detectors ' + esc(r.det.join(', ')) + '</p>';
    h += '<div class="big ' + cls + '">' + r.ps.toFixed(2) + '\\u00B0</div>';
    h += '<div class="delta ' + dcls + '">best-S/N vs official &middot; ' +
         (d === 0 ? 'unchanged from T90 window'
                  : (d < 0 ? 'improved by ' + Math.abs(d).toFixed(2) + '\\u00B0'
                           : 'worse by ' + d.toFixed(2) + '\\u00B0') + ' vs T90 window') + '</div>';
    h += '<div class="sep"></div>';
    h += '<dl class="kv">'
      + '<dt>Official centroid</dt><dd>' + r.oc[0].toFixed(2) + ', ' + r.oc[1].toFixed(2) + '</dd>'
      + '<dt>Mine, T90 window</dt><dd>' + r.ac[0].toFixed(2) + ', ' + r.ac[1].toFixed(2) + '</dd>'
      + '<dt>Mine, best S/N</dt><dd>' + r.pc[0].toFixed(2) + ', ' + r.pc[1].toFixed(2) + '</dd>'
      + '<dt>Separation, T90</dt><dd>' + r.as.toFixed(2) + '\\u00B0</dd>'
      + '<dt>Separation, best S/N</dt><dd>' + r.ps.toFixed(2) + '\\u00B0</dd>'
      + '<dt>Window, T90</dt><dd>' + r.aw.toFixed(2) + ' s</dd>'
      + '<dt>Window, best S/N</dt><dd>' + r.pw.toFixed(2) + ' s</dd>'
      + '<dt>90% area, official</dt><dd>' + r.oa + ' sq&deg;</dd>'
      + '<dt>90% area, mine</dt><dd>' + r.pa + ' sq&deg;</dd>'
      + '</dl>';
    if (r.x) {{
      var best = Math.min(r.x.do, r.x.da, r.x.dp);
      function mark(v) {{ return v === best ? ' \\u2190 closest' : ''; }}
      h += '<div class="extbox">'
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

  function show() {{
    var rec = current().filter(function (r) {{ return r.t === grbSel.value; }})[0];
    if (!rec) return;
    draw(rec);
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
