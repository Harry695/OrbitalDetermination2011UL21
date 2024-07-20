from astropy.io import fits
import numpy as np
tbl = fits.open("corr.fits")[1].data
rms_ra = 3600*(np.mean((tbl.field_ra - tbl.index_ra)**2.))**0.5
rms_dec = 3600*(np.mean((tbl.field_dec - tbl.index_dec)**2.))**0.5
print("RMS RA:", rms_ra)
print("RMS DEC:", rms_dec)