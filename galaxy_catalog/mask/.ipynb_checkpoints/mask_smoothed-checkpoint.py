from astropy.io import fits
from astropy.table import Table
from desimodel.footprint import tiles2pix  
import numpy as np
import healpy as hp

#This is Anand Raichoor's code, modified by Hernan Rincon to generate a mask for Iron

# make sure you are using the desi conda environemnt, and if needed, re-add desimodel to the python path,
# following instructions at https://desi.lbl.gov/trac/wiki/Pipeline/GettingStarted/Laptop

# read the tiles-specstatus to get the qa-validated tiles.
# tiles must contains (RA, DEC) so that tiles2pix works
# uncomment the filepath appropriate for running on NERSC or your local machine
fn = "/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv"
#fn = "tiles-specstatus.ecsv"
tiles = Table.read(fn)
#Iron end-date obtained form DESI wiki: https://desi.lbl.gov/trac/wiki/Pipeline/Releases/Iron
sel = (tiles["FAFLAVOR"] == "mainbright") & (tiles["QA"] == "good") & (tiles["LASTNIGHT"] <= 20220613)
tiles = tiles[sel] 

# read Anand's healpix map
# which stores in TILEIDS (npixels, npass) the list of tiles covering each pixel
# uncomment the filepath appropriate for running on NERSC or your local machine
#fn = "main-skymap-bright-goal.fits"
#fn = "/global/cfs/cdirs/desi/users/raichoor/main-status/skymaps/bright/main-skymap-bright-goal.fits" # 5 pass version (not used for DESIVAST V1)
fn = '/global/cfs/cdirs/desi/users/raichoor/main-status/skymaps/bright4pass/main-skymap-bright4pass-goal.fits'
hdr = fits.getheader(fn, 1)
d = fits.open(fn)[1].data # fits.open is much faster than fitsio.read...
nside = int((len(d['HPXPIXEL'])/12)**.5)
npass = d["TILEIDS"].shape[1] # npass=4 for bright
goal_ns = d["NPASS"] # number of planned tiles covering each pixel

# now count how many qa-validated tiles are done for each pixel
ns = np.nan + np.zeros(len(d))
ns[goal_ns > 0] = 0.0
for i in range(npass):
    sel = np.in1d(tiles["TILEID"], d["TILEIDS"][:, i])
    if sel.sum() > 0:
        tileids_i = tiles["TILEID"][sel]
        ipixs = tiles2pix(nside, tiles=tiles[sel])
        ns[ipixs] += 1


# fractional coverage
fracns = np.nan + np.zeros(len(d))
sel = goal_ns > 0
# For target completeness (patchy, but even target completeness)
# fracns[sel] = ns[sel] / npass 
# For survey completeness (smooth, but uneven target completeness)
fracns[sel] = ns[sel] / goal_ns[sel]

# store results in a new table
myd = Table()
for key in ["HPXPIXEL", "RA", "DEC"]:
    myd[key] = d[key]
myd["DESI"] = goal_ns > 0 # pixels covered by desi
myd["DONE"] = fracns == 1 # pixels where all tiles are obs.+qa-validated
myd["GOAL"] = goal_ns # pixels covered by desi
myd["FRAC"] = fracns # pixels where all tiles are obs.+qa-validated
myd.meta["HPXNSIDE"] = hdr["HPXNSIDE"]
myd.meta["HPXNEST"] = hdr["HPXNEST"]

# smooth the mask
# convolve with 2.8 deg sigma gaussian
smoothed_mask = hp.smoothing(myd['DONE'].astype(float), nest=True, sigma=2.8 * np.pi/180)
# choose a threshold between 0-1 for including pixels in mask
myd['DONE'][smoothed_mask < 0.338] = False
#manually remove leftover unwanted regions
theta, phi = hp.pix2ang(hp.get_nside(myd['DONE']), np.arange(len(myd['DONE'])), nest=True)
ra = np.rad2deg(phi)
dec = np.rad2deg(0.5 * np.pi - theta)
select = (ra < 40) 
myd['DONE'][select] = False

#save mask
myd.write("iron_mask_smoothed.fits")