import numpy as np
import matplotlib.pyplot as plt
from photutils import CircularAperture,aperture_photometry
from astropy.io import fits
import sys

hdu_1=fits.open("/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run07/alma_model_rotated.fits")
hdu_2=fits.open("/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run01/alma_model_rotated.fits")
data_1=hdu_1[0].data
data_2=hdu_2[0].data
data=data_1-data_2

aperture=CircularAperture((559.519,519.498),r=22)
phot_table=aperture_photometry(data, aperture)

plt.imshow(data,origin='lower')
aperture.plot(color='red',lw=1)
plt.xlim(520,600)
plt.ylim(480,560)
plt.show()

print(phot_table)
