import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import sys
import fnmatch
plt.style.use('fancy')

hdulist_ext=fits.open('ext.fits') # row: wl, col: z1,z2,z3
hdulist_abso=fits.open('abs.fits')
hdulist_sca=fits.open('sca.fits')

ext=hdulist_ext[0].data
abso=hdulist_abso[0].data
sca=hdulist_sca[0].data

Nzones=ext.shape[1]-1



############################################################
# Plotting
fig=plt.figure(figsize=(12,5))
gs=gridspec.GridSpec(1,Nzones,hspace=0.0)
for i in range(0,Nzones):
    ax=plt.subplot(gs[0,i])
    ax.plot(ext[:,0:1],ext[:,i+1:i+2],label=r"$\kappa_{\nu}^{\mathrm{ext}}$")
    ax.plot(abso[:,0:1],abso[:,i+1:i+2],label=r"$\kappa_{\nu}^{\mathrm{abs}}$")
    ax.plot(sca[:,0:1],sca[:,i+1:i+2],label=r"$\kappa_{\nu}^{\mathrm{sca}}$")
    ax.legend()
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend(fontsize=12,loc="lower left")
    ax.set_xlabel(r"$\lambda \, (\mu m)$")
    ax.set_ylim(5e-4,3e5)
    if i==0:
        ax.set_ylabel(r"Dust oppacities (cm$^2$/g(dust))")
plt.show()
fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/report:08-05-2020/nuevo.png")
    
