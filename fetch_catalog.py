"""
Fetch the Fermi GBM burst catalog (fermigbrst) from HEASARC and write a
manifest of trigger names to drive the bulk download.

Output: catalog/fermigbrst.csv  with one row per GRB.
"""
import io
import os
import urllib.parse
import urllib.request

from astropy.io.votable import parse_single_table

TAP = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
OUTDIR = "catalog"
OUTFILE = os.path.join(OUTDIR, "fermigbrst.csv")

# trigger_name is the bnYYMMDDfff string the FTP tree is keyed on.
ADQL = """
SELECT trigger_name, name, ra, dec, error_radius, trigger_time,
       t90, t90_error, t90_start, t50, t50_start,
       fluence, fluence_error, flux_64, flux_256, flux_1024,
       flnc_best_fitting_model, bcat_detector_mask, scat_detector_mask
FROM fermigbrst
ORDER BY trigger_name
"""


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    query = urllib.parse.urlencode({"QUERY": ADQL, "LANG": "ADQL"})
    url = f"{TAP}?{query}"
    print(f"Querying {TAP} ...")

    with urllib.request.urlopen(url, timeout=300) as resp:
        payload = resp.read()
    print(f"  received {len(payload) / 1e6:.1f} MB of VOTable")

    table = parse_single_table(io.BytesIO(payload)).to_table()
    print(f"  {len(table)} rows, columns: {list(table.colnames)}")

    # trigger_name comes back as bytes from the BINARY serialisation
    for col in table.colnames:
        if table[col].dtype.kind == "S":
            table[col] = [v.decode().strip() for v in table[col]]

    table.write(OUTFILE, format="ascii.csv", overwrite=True)
    print(f"Wrote {OUTFILE}")

    triggers = [t for t in table["trigger_name"] if t.startswith("bn")]
    print(f"  {len(triggers)} usable trigger names")
    print(f"  first: {triggers[0]}   last: {triggers[-1]}")

    years = sorted({t[2:4] for t in triggers})
    print(f"  years covered: 20{years[0]} - 20{years[-1]}")


if __name__ == "__main__":
    main()
