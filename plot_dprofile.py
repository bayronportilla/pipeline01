import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter
import os
import sys
import fnmatch
plt.style.use('fancy')

iprofile=np.loadtxt("/Users/users/bportilla/Documents/first_project/scripts/PDS70/surface_density/surface_density_PDS70.dat")
fprofile=np.loadtxt("../surface_density_PDS70.dat")
kprofile=np.loadtxt("/Users/users/bportilla/Documents/first_project/scripts/PDS70/keppler_dprofile.dat")

ix=np.reshape(iprofile[:,0:1],iprofile.shape[0])
iy=np.reshape(iprofile[:,1:2],iprofile.shape[0])
fx=np.reshape(fprofile[:,0:1],fprofile.shape[0])
fy=np.reshape(fprofile[:,1:2],fprofile.shape[0])
kx=np.reshape(kprofile[:,0:1],kprofile.shape[0])
ky=np.reshape(kprofile[:,1:2],kprofile.shape[0])

fsize=14
fig=plt.figure()
ax=plt.axes()
ax.plot(ix,iy,color="lightblue",linewidth=2.0)
ax.plot(fx,fy,color="salmon",linewidth=2.0)
ax.plot(kx,ky,color="seagreen",linewidth=2.0)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"Radial distance (AU)",fontsize=fsize)
ax.set_ylabel(r"$\Sigma_{\mathrm{dust}}$ (g/cm^2)",fontsize=fsize)
rotn=-20.0
ax.annotate("initial surface density",xy=(0.1,6.9),ha='left',va='top',rotation=rotn,color="grey")
ax.annotate("modified surface density",xy=(0.05,0.07),ha='left',va='top',rotation=rotn,color="grey")
ax.annotate("Keppler et al. 2018",xy=(0.06,0.24),ha='left',va='top',rotation=rotn,color="grey")
ax.axvline(40.0,1e-5,10,linestyle="--",color="lightgrey")
ax.annotate(r"$R_{\mathrm{tap}}$",xy=(43.0,1.7),ha='left',va='top',color="grey")
#ax.minorticks_off()
ax.tick_params(labelsize=14)
ax.set_xlim(min(fx),max(fx))
ax.set_ylim(1e-5,)


#plt.show()
#sys.exit()
name="surface_density"
fig.savefig("../%s.png"%(name))
fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/paper_figures/%s.png"%(name))


