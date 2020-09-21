import numpy as np
import matplotlib.pyplot as plt
from gofish import imagecube
import sys
import matplotlib.gridspec as gridspec


############################################################
# Load data
cube=imagecube('../alma_model_rotated.fits')


############################################################
# Main routine
def radial_cut(PA_disk,inc,d,padir,widir,dr,x_off,y_off):
    PAmin=padir-0.5*widir
    PAmax=padir+0.5*widir
    x,y,dy=cube.radial_profile(inc=inc, PA=PA_disk, 
                               dist=d, 
                               PA_min=PAmin,PA_max=PAmax,dr=dr,
                               assume_correlated=False,x0=x_off,y0=y_off)
    return x,y,dy


############################################################
# Input parameters
ddisk=113.43 # Distance (pc)
drsample=0.02 # Step size (arcsec)
padisk=158.6 # Disk position angle (deg)
incdisk=49.7 # Disk inclination (deg)
padir=118.4 # Direction of interest, w.r.t. red shifted semi-major axis (deg)
widir=20.0 # width of the cone (deg)
nbins=1 # Number of cones (odd integer)


############################################################
# Plotting
fig=plt.figure()
for i in range(int(-0.5*(nbins-1)),int(0.5*(nbins-1))+1):
    angle=padir+i*widir
    x,y,dy=radial_cut(padisk,incdisk,ddisk,angle,widir,drsample,0.0,0.0)
    plt.errorbar(x*ddisk,y,dy,fmt=".",markersize=4,elinewidth=1,
                 label="PA=%.1f deg"%(angle+padisk))
    
plt.axvspan(18.0,24.0,color='red',alpha=0.1)
plt.xlabel("r (AU)")
plt.ylabel("Flux density (mJy/beam)")
plt.legend()
plt.xlim(0,120)
plt.show()
#fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/reports/report:11-09-2020/coneobs_b.png")

############################################################
# Writing file
f=open("../cut_along_c_mod.dat","w")
for (i,j,k) in zip(x,y,dy):
    f.write("%.15e %.15e %.15e\n"%(i*ddisk,j,k)) # AU,mJy/beam,mJy/beam
f.close()
sys.exit()




