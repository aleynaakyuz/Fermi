# Fermi GBM burst localization

Reproducing the Fermi Gamma-ray Burst Monitor's own gamma-ray burst positions
from archival data, and checking the result against both GBM's published sky
maps and independent positions from other observatories.

**[Browse the results →](https://aleynaakyuz.github.io/Fermi/)**

## What this does

GBM locates a burst by comparing counts across its twelve NaI detectors, whose
responses point in different directions. This project re-runs that calculation
independently for every GBM burst that has an official sky map — 2,049 of them
— twice each:

- using the catalog's full **T90** interval as the source window, and
- using only the **brightest sub-interval**, the way a human Burst Advocate
  would pick it.

Each result is compared against the official Fermi HEALPix map and, for the 284
bursts where Swift or the Fermi LAT supplied an arcsecond position, against
something independent of GBM entirely.

## Main results

| | Long (1,695) | Short (354) |
|---|---|---|
| Median separation from official, full T90 | 4.77° | 6.48° |
| Median separation, best signal-to-noise window | **3.20°** | **4.70°** |
| Bursts improved by the narrower window | 1,078 (64%) | 185 (52%) |

Trimming the source window to the bright part of the burst helps both classes
(Wilcoxon signed-rank, p < 10⁻³ for short and p < 10⁻⁵⁷ for long). It helps
long bursts more reliably, which fits the mechanism: trimming discards time
carrying background but little signal, and a short burst has hardly any such
time to discard.

Measured against the 284 bursts with an independent position, the reproduction
is **as accurate as the official pipeline** — closer to the truth on 48% of
them, a coin flip. Separation from the official map measures agreement with
GBM, not accuracy; the official map is itself a median 4.31° from truth for
those bursts.

## Pages

| Page | Contents |
|---|---|
| [Burst viewer](https://aleynaakyuz.github.io/Fermi/viewer.html) | Pick a class and a GRB; see all four localizations plotted together |
| [Long GRBs](https://aleynaakyuz.github.io/Fermi/long.html) | 1,695 bursts, summary charts, full sortable table |
| [Short GRBs](https://aleynaakyuz.github.io/Fermi/short.html) | 354 bursts, including whether trimming helps them |
| [Cross-mission](https://aleynaakyuz.github.io/Fermi/external.html) | The 969 GBM bursts with a position from another instrument |
| [Replication report](https://aleynaakyuz.github.io/Fermi/replication.html) | The original ten-burst study and what closed the gap |

## Reproducing

Requires `astro-gdt-fermi` and `astropy`.

```bash
python fetch_catalog.py          # burst catalog from HEASARC (4,400 GRBs)
python download_gbm.py           # trigdat + official sky maps (5.2 GB)
python external_localizations.py # cross-match Swift / LAT catalogs
python run_all.py                # localize everything, both windows
python extract_contours.py       # contour polygons for the viewer
python build_class_pages.py      # long.html, short.html
python build_viewer.py           # viewer.html
python build_external_page.py    # external.html
python publish_github.py         # wrap pages into docs/ for GitHub Pages
```

`download_gbm.py` and the analysis scripts are resumable — re-running skips
work already done, so an interrupted run picks up where it stopped.

The 5.2 GB of burst data is **not** in this repository; `download_gbm.py`
fetches it from the HEASARC archive.

## Scripts

| Script | Purpose |
|---|---|
| `fetch_catalog.py` | HEASARC TAP query for the GBM burst catalog |
| `download_gbm.py` | Resumable bulk download of trigdat and HEALPix maps |
| `external_localizations.py` | Cross-match against Swift and Fermi LAT catalogs |
| `analyze_grbs.py` | The localization itself: windowing, background, DoL fit |
| `run_all.py` | Parallel run over the full catalog, both window modes |
| `extract_contours.py` | 90% contour polygons, so the viewer can draw any burst |
| `build_*.py` | Page generation |
| `sweep_erange.py` | Energy-band parameter sweep |

## Notes on method

- **Energy band** is the standard 50–300 keV. A sweep over all eight bands
  trigdat's channel structure allows found narrowing makes things markedly
  worse: restricting to a single channel roughly triples the error.
- **`legacy_DoL` expects count rates, not counts**, despite its docstring. The
  distinction is harmless for a short window, where the error cancels in the
  flux normalisation, but over a 100 s window it inflates every chi-square past
  the routine's internal clamp and the fit fails outright.
- **Separations are quantized** at roughly 0.46°, the HEALPix NSIDE 128 pixel
  scale, because a map's centroid is reported as its highest-probability pixel
  centre rather than a weighted mean.
- **The systematic model cannot be matched per burst.** Every official HEALPix
  file carries `LMETHOD = 'Interactive'`, written unconditionally, so the
  header does not record whether a human or RoboBA produced the map. All
  bursts here use `hitl_model`.
- **12 of 2,061 bursts failed**: 7 where the grid search found no acceptable
  spectrum, 4 where the data span left no room for a clean background window,
  and 1 singular matrix in the background fit.

## Data

All inputs are public NASA archival data: the Fermi GBM burst catalog
(`fermigbrst`), the Fermi LAT GRB catalog (`fermilgrb`), and the Swift GRB
catalogs (`swiftgrb`, `swiftgrbba`), all from
[HEASARC](https://heasarc.gsfc.nasa.gov/), plus GBM trigger data from the
Fermi archive.
