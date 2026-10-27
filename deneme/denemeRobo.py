"""
Recreate the Fermi GBM ground localization of a known short GRB
(GRB 200826A, trigger bn200826187) using the DoL code in gdt-fermi,
then compare it with the official GBM map and the optical afterglow.

Follows the DoL tutorial in the gdt-fermi documentation:
https://astro-gdt.readthedocs.io/projects/astro-gdt-fermi/en/latest/missions/fermi/gbm/localization/dol/dol.html

Setup:
    pip install astro-gdt-fermi
    gdt-data init

Workflow:
    1. Run once. Open lightcurve.png and check that SRC_TIME covers the
       burst and BG_TIMES sit on flat background either side.
    2. Adjust the windows below if needed and run again.
"""
import glob
import os
import urllib.request

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

from gdt.missions.fermi.gbm.finders import TriggerFinder
from gdt.missions.fermi.gbm.trigdat import Trigdat
from gdt.missions.fermi.time import Time
from gdt.core.plot.lightcurve import Lightcurve
from gdt.core.background.fitter import BackgroundFitter
from gdt.core.background.binned import Polynomial
from gdt.core.plot.sky import EquatorialPlot
from gdt.missions.fermi.gbm.localization import GbmHealPix, hitl_model, robo_ba_model
from gdt.missions.fermi.gbm.localization.dol.legacy_dol import legacy_DoL

# ---------------------------------------------------------------- settings
BURST = "200826187"                    # GBM burst number (bnYYMMDDfff)
SRC_TIME = (-0.256, 1.792)             # source window, s from trigger (64 ms multiples)
BG_TIMES = [(-8.0, -1.5), (4.096, 14.0)]  # background windows, s from trigger
LOC_ERANGE = (50.0, 300.0)             # standard GBM localization band, keV
TIMESCALE = 64                         # trigdat binning in ms (short GRB)
GRB_TYPE = "short"                     # RoboBA systematic model: "short" or "long"

# Reference positions to compare against
TRUTH = SkyCoord(6.7855741, 34.0273043, unit="deg")    # ZTF20abwysqy afterglow (GCN 28295)
ROBOBA = SkyCoord(4.7, 35.3, unit="deg")               # GBM automated position (GCN 28284)

OFFICIAL_MAP_URL = ("https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2020/"
                    "bn200826187/quicklook/glg_healpix_all_bn200826187.fit")


def report(label, healpix):
    """Print centroid, 90% area and where the afterglow falls in a map."""
    ra, dec = healpix.centroid
    centroid = SkyCoord(ra, dec, unit="deg")
    print(f"  {label}")
    print(f"    centroid           RA {ra:7.2f}  Dec {dec:+6.2f}")
    print(f"    90% area           {healpix.area(0.9):8.1f} sq deg")
    print(f"    afterglow offset   {centroid.separation(TRUTH).deg:6.2f} deg")
    print(f"    afterglow lies at  {100 * healpix.confidence(TRUTH.ra.deg, TRUTH.dec.deg):5.1f}% contour")


# ------------------------------------------------------- 1. get trigdat data
TriggerFinder(BURST).get_trigdat(".")
trigdat_file = sorted(glob.glob(f"glg_trigdat_all_bn{BURST}_v*.fit"))[-1]
trigdat = Trigdat.open(trigdat_file)
print(f"Opened {trigdat_file}")
print(f"Triggered detectors: {trigdat.triggered_detectors}")

# --------------------------------------------- 2. source and background counts
src_counts, src_exposure = [], []
bg_counts, bg_exposure = [], []
bg_trigdet = []

for det in trigdat._detectors:
    phaii = trigdat.to_phaii(det, timescale=TIMESCALE)

    src_bin = phaii.data.integrate_time(*SRC_TIME)
    src_counts.append(src_bin.counts.astype(np.int32))
    src_exposure.append(src_bin.exposure)

    fitter = BackgroundFitter.from_phaii(phaii, Polynomial, time_ranges=BG_TIMES)
    fitter.fit(order=1)
    bg_rates = fitter.interpolate_bins(phaii.data.tstart, phaii.data.tstop)
    bg_bin = bg_rates.integrate_time(*SRC_TIME)
    bg_counts.append(bg_bin.counts.astype(np.int32))
    bg_exposure.append(bg_bin.exposure)

    if det in trigdat.triggered_detectors:
        bg_trigdet.append(bg_rates.slice_energy(*LOC_ERANGE))

    net = int(src_counts[-1].sum() - bg_counts[-1].sum())
    print(f"  {det}: net counts in source window = {net}")

avg_src_exposure = np.mean(src_exposure)
avg_bg_exposure = np.mean(bg_exposure)

# Light curve of the triggered detectors with the windows marked
summed = trigdat.sum_detectors(trigdat.triggered_detectors, timescale=TIMESCALE)
summed = summed.slice_energy(LOC_ERANGE).slice_time((-8, 14))
lc = Lightcurve(summed.to_lightcurve())
lc.set_background(bg_trigdet[0].sum_time(bg_trigdet))
ax = plt.gca()
ax.axvspan(*SRC_TIME, color="C1", alpha=0.2)
for window in BG_TIMES:
    ax.axvspan(*window, color="C0", alpha=0.2)
plt.savefig("lightcurve.png", dpi=150, bbox_inches="tight")
plt.close("all")
print("Saved lightcurve.png (orange = source, blue = background)")

# ------------------------------------------ 3. energy bins and spacecraft state
phaii = trigdat.to_phaii(trigdat._detectors[0], timescale=TIMESCALE)
energies = np.concatenate([phaii.data.emin, [phaii.data.emax[-1]]])
crange = [np.digitize(e, energies, right=True) - 1 for e in LOC_ERANGE]

tcenter = trigdat.trigtime + 0.5 * sum(SRC_TIME)
frame = trigdat.poshist.at(Time(tcenter, format="fermi"))
scpos = frame.obsgeoloc.xyz.to_value("km")
quaternion = frame.quaternion.scalar_last

flight_ra = trigdat.headers["PRIMARY"]["RA_OBJ"]
flight_dec = trigdat.headers["PRIMARY"]["DEC_OBJ"]
flight = SkyCoord(flight_ra, flight_dec, unit="deg")

# --------------------------------------------------------- 4. run the DoL fit
dol = legacy_DoL()
loc = dol.eval(crange, np.array(src_counts), np.array(bg_counts),
               avg_src_exposure, avg_bg_exposure,
               scpos, quaternion, energies, flight_ra, flight_dec, int(tcenter),
               scat_opt=1)

best_ra, best_dec = np.degrees([loc["best"]["ra"], loc["best"]["dec"]])
best = SkyCoord(best_ra, best_dec, unit="deg")
az = np.degrees(loc["best"]["az"])

print("\nBest-fit positions")
print(f"  Your DoL fit       RA {best_ra:7.2f}  Dec {best_dec:+6.2f}"
      f"   (68% stat radius {loc['best']['err']:.2f} deg)")
print(f"  Flight software    RA {flight_ra:7.2f}  Dec {flight_dec:+6.2f}")
print(f"  GBM RoboBA (GCN)   RA {ROBOBA.ra.deg:7.2f}  Dec {ROBOBA.dec.deg:+6.2f}")
print(f"  Afterglow (truth)  RA {TRUTH.ra.deg:7.2f}  Dec {TRUTH.dec.deg:+6.2f}")
print("\nOffsets from the afterglow")
print(f"  Your DoL fit       {best.separation(TRUTH).deg:6.2f} deg")
print(f"  Flight software    {flight.separation(TRUTH).deg:6.2f} deg")
print(f"  GBM RoboBA         {ROBOBA.separation(TRUTH).deg:6.2f} deg")

# -------------------------------------- 5. probability maps with systematics
# The same statistical map, convolved with two different systematic models:
#   HitL:   depends on the burst's spacecraft azimuth (human Burst Advocate)
#   RoboBA: two-Gaussian mixture from Goldstein et al. 2020, ApJ 895, 40,
#           with separate parameters for short and long GRBs
# GbmHealPix.convolve passes the spacecraft quaternion and position itself;
# the extra arguments go to the model function.
stat_map = dol.to_GbmHealPix(loc, frame)
hitl_map = stat_map.convolve(hitl_model, az).remove_earth()
roboba_map = stat_map.convolve(robo_ba_model, GRB_TYPE).remove_earth()

print("\nMap comparison")
report("Your map, statistical only", stat_map)
report("Your map, HitL systematics", hitl_map)
report(f"Your map, RoboBA systematics ({GRB_TYPE})", roboba_map)

official_file = os.path.basename(OFFICIAL_MAP_URL)
try:
    if not os.path.exists(official_file):
        urllib.request.urlretrieve(OFFICIAL_MAP_URL, official_file)
    official = GbmHealPix.open(official_file).remove_earth()
    report("Official GBM map (RoboBA)", official)
except Exception as err:
    official = None
    print(f"  Could not load the official map ({err}); download it from\n  {OFFICIAL_MAP_URL}")

# --------------------------------------------------------------- 6. sky plots
maps = [("skymap_hitl.png", hitl_map), ("skymap_roboba.png", roboba_map)]
if official is not None:
    maps.append(("skymap_official.png", official))
for filename, healpix in maps:
    plot = EquatorialPlot()
    plot.add_localization(healpix, clevels=[0.90, 0.50], gradient=False)
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close("all")
    print(f"Saved {filename}")
