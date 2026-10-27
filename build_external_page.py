"""
Render catalog/external_localizations.csv as a browsable, filterable page.

969 rows is not readable as a flat dump, so the page is built as a tool:
filter by instrument, search by name, sort any column. Pure ASCII output,
self-contained, no external assets.
"""
import csv
import collections
import datetime
import html
import json
import os

SRC = "catalog/external_localizations.csv"
OUT = "site/external.html"

OFFICIAL = "#2f6fd0"


def clean(value):
    value = (value or "").strip()
    return "" if value in ("--", "nan", "None", "NULL") else value


def main():
    rows = list(csv.DictReader(open(SRC, newline="")))

    records, counts = [], collections.Counter()
    years = collections.Counter()
    obs_counts = collections.Counter()
    n_redshift = 0

    for r in rows:
        attribution = r["attribution"]
        counts[attribution] += 1
        year = "20" + r["trigger_name"][2:4]
        years[year] += 1
        obs = clean(r.get("swift_trigger_obs"))
        if obs:
            obs_counts[obs] += 1
        redshift = clean(r.get("swift_redshift"))
        if redshift:
            n_redshift += 1

        if r["swift_name"]:
            ext_name = clean(r["swift_name"])
            ext_ra, ext_dec = clean(r["swift_ra"]), clean(r["swift_dec"])
        elif r["lat_name"]:
            ext_name = clean(r["lat_name"])
            ext_ra, ext_dec = clean(r["lat_ra"]), clean(r["lat_dec"])
        else:
            ext_name, ext_ra, ext_dec = "", "", ""

        records.append({
            "t": r["trigger_name"],
            "y": year,
            "t90": clean(r["t90"]),
            "ra": clean(r["catalog_ra"]),
            "dec": clean(r["catalog_dec"]),
            "err": clean(r["catalog_err_deg"]),
            "a": attribution,
            "en": ext_name,
            "era": ext_ra,
            "edec": ext_dec,
            "off": clean(r["catalog_vs_external_deg"]),
            "z": redshift,
            "obs": obs,
            "q": r.get("match_quality", "ok"),
        })

    total = len(records)
    n_suspect = sum(1 for r in records if r["q"] == "suspect")
    n_own = sum(1 for r in records if r["q"] == "gbm_own_position")
    swift_any = counts["Swift"] + counts["Swift+Fermi-LAT"]
    lat_any = counts["Fermi-LAT"] + counts["Swift+Fermi-LAT"]
    built = datetime.date.today().isoformat()
    data = json.dumps(records, separators=(",", ":"))

    year_rows = "".join(
        f"<tr><th scope=\"row\">{y}</th><td class=\"num\">{n}</td></tr>"
        for y, n in sorted(years.items()))
    obs_rows = "".join(
        f"<tr><th scope=\"row\">{html.escape(k)}</th><td class=\"num\">{v}</td></tr>"
        for k, v in obs_counts.most_common())

    page = f"""<title>Externally Localized GBM Bursts</title>
<style>
  :root {{
    --paper:#f5f7fa; --surface:#ffffff; --surface-2:#eef1f6;
    --ink:#141b26; --ink-soft:#46525f; --muted:#6b7686; --rule:#d8dee8;
    --accent:{OFFICIAL}; --swift:#2f6fd0; --lat:#9a5bc4; --none:#8a7a3f;
    --serif: ui-serif, Georgia, "Iowan Old Style", "Times New Roman", serif;
    --sans: system-ui, -apple-system, "Segoe UI", Roboto, sans-serif;
    --mono: ui-monospace, SFMono-Regular, "SF Mono", Menlo, Consolas, monospace;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
      --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
      --accent:#6fa3ea; --swift:#6fa3ea; --lat:#bb8ade; --none:#c6ad63;
    }}
  }}
  :root[data-theme="dark"] {{
    --paper:#0e131b; --surface:#151c27; --surface-2:#1b2430;
    --ink:#e8edf4; --ink-soft:#b3becd; --muted:#8995a6; --rule:#2a3545;
    --accent:#6fa3ea; --swift:#6fa3ea; --lat:#bb8ade; --none:#c6ad63;
  }}
  *,*::before,*::after {{ box-sizing:border-box; }}
  body {{
    margin:0; background:var(--paper); color:var(--ink);
    font-family:var(--sans); font-size:16px; line-height:1.6;
    -webkit-font-smoothing:antialiased;
  }}
  .wrap {{
    max-width:1140px; margin:0 auto;
    padding:clamp(1.5rem,4vw,3rem) clamp(.8rem,3vw,1.6rem) 4rem;
    display:flex; flex-direction:column; gap:2.2rem;
  }}
  .eyebrow {{
    font-family:var(--mono); font-size:.72rem; letter-spacing:.14em;
    text-transform:uppercase; color:var(--muted); margin:0;
  }}
  h1 {{
    font-family:var(--serif); font-weight:600; margin:.5rem 0 0;
    font-size:clamp(1.7rem,4vw,2.5rem); line-height:1.14;
    letter-spacing:-.012em; text-wrap:balance;
  }}
  .standfirst {{
    font-family:var(--serif); color:var(--ink-soft); max-width:64ch;
    font-size:clamp(1rem,2vw,1.12rem); margin:.9rem 0 0;
  }}
  h2 {{
    font-family:var(--serif); font-size:1.3rem; font-weight:600;
    margin:0 0 .8rem; text-wrap:balance;
  }}
  p {{ max-width:68ch; margin:0 0 .9rem; }}
  code {{ font-family:var(--mono); font-size:.9em; }}

  .tiles {{
    display:grid; grid-template-columns:repeat(auto-fit,minmax(170px,1fr));
    gap:1px; background:var(--rule); border:1px solid var(--rule);
    border-radius:10px; overflow:hidden;
  }}
  .tile {{ background:var(--surface); padding:1rem 1.15rem; }}
  .tile .v {{
    font-family:var(--serif); font-size:1.8rem; line-height:1.1; display:block;
    font-variant-numeric:tabular-nums;
  }}
  .tile .k {{ font-size:.78rem; color:var(--muted); display:block; margin-top:.3rem; }}

  .controls {{
    display:flex; flex-wrap:wrap; gap:.6rem; align-items:center;
    padding:.9rem 1rem; background:var(--surface);
    border:1px solid var(--rule); border-radius:10px;
  }}
  .controls input[type="search"] {{
    flex:1 1 220px; min-width:170px; padding:.5rem .7rem;
    border:1px solid var(--rule); border-radius:7px;
    background:var(--paper); color:var(--ink);
    font-family:var(--sans); font-size:.9rem;
  }}
  .controls input:focus-visible, .pill:focus-visible {{
    outline:2px solid var(--accent); outline-offset:2px;
  }}
  .pill {{
    font:inherit; font-size:.82rem; cursor:pointer;
    padding:.4rem .75rem; border-radius:999px;
    border:1px solid var(--rule); background:var(--paper); color:var(--ink-soft);
  }}
  .pill[aria-pressed="true"] {{
    background:var(--accent); border-color:transparent; color:#fff; font-weight:600;
  }}
  .count {{ font-size:.82rem; color:var(--muted); margin-left:auto; }}

  .table-scroll {{
    overflow-x:auto; border:1px solid var(--rule); border-radius:10px;
    max-height:70vh; overflow-y:auto;
  }}
  table {{ border-collapse:collapse; width:100%; min-width:900px;
           background:var(--surface); font-size:.86rem; }}
  th,td {{ padding:.42rem .7rem; border-bottom:1px solid var(--rule); text-align:left; }}
  thead th {{
    position:sticky; top:0; z-index:2; background:var(--surface-2);
    font-size:.7rem; letter-spacing:.06em; text-transform:uppercase;
    color:var(--muted); font-weight:650; white-space:nowrap; cursor:pointer;
    user-select:none;
  }}
  thead th::after {{ content:" \\2195"; opacity:.3; }}
  thead th.asc::after {{ content:" \\2191"; opacity:.9; }}
  thead th.desc::after {{ content:" \\2193"; opacity:.9; }}
  td.num {{ text-align:right; font-variant-numeric:tabular-nums; }}
  tbody td:first-child {{ font-family:var(--mono); font-size:.82rem; white-space:nowrap; }}
  tbody tr:nth-child(even) {{ background:color-mix(in srgb, var(--surface-2) 45%, transparent); }}

  .src {{ font-weight:600; white-space:nowrap; }}
  .src.Swift {{ color:var(--swift); }}
  .src.LAT {{ color:var(--lat); }}
  .src.Both {{ color:var(--lat); }}
  .src.None {{ color:var(--none); }}
  td.suspect {{ color:#bf3b36; font-weight:650; }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) td.suspect {{ color:#e8706a; }}
  }}
  :root[data-theme="dark"] td.suspect {{ color:#e8706a; }}

  .side {{
    display:grid; grid-template-columns:repeat(auto-fit,minmax(260px,1fr)); gap:1.4rem;
  }}
  .side table {{ min-width:0; font-size:.82rem; }}
  .side .table-scroll {{ max-height:none; }}
  .callout {{
    background:var(--surface-2); border-left:3px solid var(--accent);
    border-radius:0 8px 8px 0; padding:.9rem 1.1rem; margin:0 0 1rem;
  }}
  .callout p:last-child {{ margin-bottom:0; }}
  footer {{
    border-top:1px solid var(--rule); padding-top:1.2rem;
    color:var(--muted); font-size:.82rem;
  }}
</style>

<div class="wrap">

  <header>
    <p class="eyebrow">Fermi GBM &middot; cross-mission localizations</p>
    <h1>GBM bursts that somebody else pinned down</h1>
    <p class="standfirst">
      Of the 4,400 GRBs in the Fermi GBM burst catalog, {total} carry a position
      that GBM could not have produced on its own. GBM's ground localization
      bottoms out at a systematic floor of 2&ndash;4 degrees, so every burst
      here was located by another instrument &mdash; usually Swift, sometimes
      the Fermi LAT, occasionally something that only shows up in a circular.
    </p>
  </header>

  <div class="tiles">
    <div class="tile"><span class="v">{total}</span><span class="k">bursts with an external position</span></div>
    <div class="tile"><span class="v">{swift_any}</span><span class="k">matched to Swift</span></div>
    <div class="tile"><span class="v">{lat_any}</span><span class="k">matched to Fermi-LAT</span></div>
    <div class="tile"><span class="v">{n_redshift}</span><span class="k">with a measured redshift</span></div>
  </div>

  <section>
    <h2>The list</h2>
    <div class="controls">
      <input type="search" id="q" placeholder="Search trigger, GRB name, year..."
             aria-label="Search the list" />
      <button class="pill" data-f="all" aria-pressed="true">All {total}</button>
      <button class="pill" data-f="Swift" aria-pressed="false">Swift {counts['Swift']}</button>
      <button class="pill" data-f="Fermi-LAT" aria-pressed="false">LAT {counts['Fermi-LAT']}</button>
      <button class="pill" data-f="Swift+Fermi-LAT" aria-pressed="false">Both {counts['Swift+Fermi-LAT']}</button>
      <button class="pill" data-f="unattributed" aria-pressed="false">Unattributed {counts['unattributed']}</button>
      <button class="pill" data-f="z" aria-pressed="false">Has redshift</button>
      <button class="pill" data-f="suspect" aria-pressed="false">Suspect {n_suspect}</button>
      <span class="count" id="count"></span>
    </div>

    <div class="table-scroll" style="margin-top:.9rem">
      <table id="tbl">
        <thead>
          <tr>
            <th data-k="t">Trigger</th>
            <th data-k="t90" class="num">T90 (s)</th>
            <th data-k="ra" class="num">RA</th>
            <th data-k="dec" class="num">Dec</th>
            <th data-k="err" class="num">Err (deg)</th>
            <th data-k="a">Localized by</th>
            <th data-k="en">External name</th>
            <th data-k="obs">Triggered by</th>
            <th data-k="off" class="num">Offset (deg)</th>
            <th data-k="z" class="num">z</th>
          </tr>
        </thead>
        <tbody id="body"></tbody>
      </table>
    </div>
  </section>

  <section>
    <h2>How these were identified</h2>
    <div class="callout">
      <p>
        <strong>Flagged</strong> by error radius: the GBM catalog adopts an
        external position when one exists, and GBM alone cannot reach better
        than a couple of degrees. Anything listed below 0.5&deg; is carrying
        somebody else's answer. <strong>Attributed</strong> by cross-match: the
        Fermi LAT catalog joins directly on its GBM designation; Swift is
        matched on trigger time within 60&thinsp;s.
      </p>
    </div>
    <p>
      The two signals agree well. Swift matches sit a median 0.0056&deg;
      (20 arcsec) from the position the GBM catalog lists, with 92% inside
      0.1&deg; &mdash; which is what it looks like when the catalog has simply
      adopted the Swift position outright.
    </p>
    <p>
      LAT is different: only 36% of LAT matches land within 0.1&deg;, because
      the GBM catalog generally does <em>not</em> adopt LAT positions. Those
      rows are genuine joint detections where the listed position is still
      GBM's own. Worth keeping in mind if you use this as a truth set &mdash;
      a LAT match does not guarantee a precise catalog position.
    </p>
    <p>
      Six matches are flagged <strong>suspect</strong>: the catalog carries a
      precise adopted position, but the burst found by cross-match sits more
      than a degree away, which means the match caught the wrong burst rather
      than that the position is wrong. <code>bn160625945</code> is the clearest
      case &mdash; GRB 160625B paired with Swift's separate GRB 160625A. A
      further {n_own} rows show large offsets for an innocent reason: the
      catalog holds GBM's own position there, so the offset is simply GBM's
      localization error.
    </p>
    <p>
      The {counts['unattributed']} unattributed rows are all error-radius
      flagged, so they do have external positions, from instruments not queried
      here or from afterglows reported only in GCN circulars. They cluster
      heavily in recent years, where the cross-match catalogs have not caught
      up yet.
    </p>
  </section>

  <section>
    <h2>Breakdown</h2>
    <div class="side">
      <div>
        <h3 style="font-size:.78rem;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:0 0 .5rem">By year</h3>
        <div class="table-scroll">
          <table><tbody>{year_rows}</tbody></table>
        </div>
      </div>
      <div>
        <h3 style="font-size:.78rem;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:0 0 .5rem">Which observatory triggered (Swift-matched rows)</h3>
        <div class="table-scroll">
          <table><tbody>{obs_rows}</tbody></table>
        </div>
      </div>
    </div>
  </section>

  <footer>
    Built {built} from the HEASARC <code>fermigbrst</code>, <code>fermilgrb</code>,
    <code>swiftgrb</code> and <code>swiftgrbba</code> catalogs.
    Full data in <code>catalog/external_localizations.csv</code>.
  </footer>
</div>

<script>
(function () {{
  var DATA = {data};
  var body = document.getElementById('body');
  var countEl = document.getElementById('count');
  var q = document.getElementById('q');
  var filter = 'all', sortKey = 't', sortDir = 1;

  var NUM = {{t90:1, ra:1, dec:1, err:1, off:1, z:1}};

  function label(a) {{
    if (a === 'Swift') return ['Swift', 'Swift'];
    if (a === 'Fermi-LAT') return ['Fermi-LAT', 'LAT'];
    if (a === 'Swift+Fermi-LAT') return ['Swift + LAT', 'Both'];
    return ['unattributed', 'None'];
  }}

  function esc(s) {{
    return String(s == null ? '' : s).replace(/[&<>"]/g, function (c) {{
      return {{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}}[c];
    }});
  }}

  function visible() {{
    var term = q.value.trim().toLowerCase();
    return DATA.filter(function (r) {{
      if (filter === 'z') {{ if (!r.z) return false; }}
      else if (filter === 'suspect') {{ if (r.q !== 'suspect') return false; }}
      else if (filter !== 'all' && r.a !== filter) return false;
      if (!term) return true;
      return (r.t + ' ' + r.en + ' ' + r.y + ' ' + r.obs).toLowerCase().indexOf(term) >= 0;
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

    var html = rows.map(function (r) {{
      var L = label(r.a);
      return '<tr>'
        + '<td>' + esc(r.t) + '</td>'
        + '<td class="num">' + (r.t90 ? parseFloat(r.t90).toFixed(2) : '') + '</td>'
        + '<td class="num">' + (r.ra ? parseFloat(r.ra).toFixed(3) : '') + '</td>'
        + '<td class="num">' + (r.dec ? parseFloat(r.dec).toFixed(3) : '') + '</td>'
        + '<td class="num">' + (r.err !== '' ? parseFloat(r.err).toFixed(3) : '') + '</td>'
        + '<td class="src ' + L[1] + '">' + L[0] + '</td>'
        + '<td>' + esc(r.en) + '</td>'
        + '<td>' + esc(r.obs) + '</td>'
        + '<td class="num' + (r.q === 'suspect' ? ' suspect' : '') + '">'
        + (r.off ? parseFloat(r.off).toFixed(3) : '')
        + (r.q === 'suspect' ? ' !' : '') + '</td>'
        + '<td class="num">' + esc(r.z) + '</td>'
        + '</tr>';
    }}).join('');

    body.innerHTML = html;
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

    os.makedirs("site", exist_ok=True)
    with open(OUT, "w") as handle:
        handle.write(page)
    non_ascii = sum(1 for ch in page if ord(ch) > 127)
    print(f"wrote {OUT} ({os.path.getsize(OUT)/1e6:.2f} MB, {total} rows, "
          f"non-ascii={non_ascii})")


if __name__ == "__main__":
    main()
