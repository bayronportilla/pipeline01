import numpy as np
import matplotlib.pyplot as plt
from gofish import imagecube
import sys
import matplotlib.gridspec as gridspec
from astropy.io import fits
from astropy.coordinates import Angle
import astropy.units as u
"""
cube=imagecube('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',FOV=3)
#cube2=imagecube('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',FOV=3)
#cube2=imagecube('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70.fits',FOV=3)
cube2=imagecube('/data/users/bportilla/runs/final_runs/1BS/run04_copy/alma_model_rotated.fits',FOV=3)
xm, ym, dym = cube.radial_profile(inc=49.7,PA=158.6,dist=113.43,x0=0.02,y0=0.02)
xm2, ym2, dym2 = cube2.radial_profile(inc=49.7,PA=158.6,dist=113.43)
#x=[i*]
#print(xm2)
#sys.exit()
plt.errorbar(xm,ym*1000,dym*1000,fmt='.',label="obs")
plt.errorbar(xm2,ym2,dym2,fmt='.',label="mod")
plt.legend()
plt.show()
sys.exit()
f=open("gofishdat.dat","w")
for (i,j,k) in zip(xm,ym,dym):
    f.write("%.13e %.13e %.13e \n"%(i*113.43,j/np.nanmax(ym),k/np.nanmax(ym)))
f.close()
plt.errorbar(xm*113.43,ym/np.nanmax(ym),dym/np.nanmax(ym))
plt.show()
"""

def get_profile(file,PA_disk,inc,d,padir,widir,dr):
    cube=imagecube(file,FOV=3)
    PAmin=padir-0.5*widir
    PAmax=padir+0.5*widir
    x,y,dy=cube.radial_profile(inc=inc, PA=PA_disk, 
                               dist=d, 
                               PA_min=PAmin,PA_max=PAmax,dr=dr,
                               assume_correlated=False)
                               #exclude_PA=False)
    return x,y,dy

ddisk=113.43
drsample=0.02
padisk=158.6
incdisk=49.7
padir=118.4
#padir=-61.6
widir=20.0
nbins=1 # odd integer

fig=plt.figure()
for i in range(int(-0.5*(nbins-1)),int(0.5*(nbins-1))+1):
    angle=padir+i*widir
    x,y,dy=get_profile('/data/users/bportilla/runs/final_runs/1BS/run01_copy/alma_model_rotated.fits',
                       padisk,incdisk,ddisk,angle,widir,drsample)
    x2,y2,dy2=get_profile('/data/users/bportilla/runs/final_runs/1BS/run04_copy/alma_model_rotated.fits',
                       padisk,incdisk,ddisk,angle,widir,drsample)
    """
    plt.errorbar(x*ddisk,y/np.nanmax(ym),dy/np.nanmax(ym),fmt=".",markersize=4,elinewidth=1,
                 label="PA=%.1f deg"%(angle+padisk))
    """
    plt.errorbar(x*ddisk,y,dy,fmt=".",markersize=4,elinewidth=1,
                 label="no planet")
    plt.errorbar(x2*ddisk,y2,dy2,fmt=".",markersize=4,elinewidth=1,
                 label="with planet")
    

plt.axvspan(34.1,43.1,color='red',alpha=0.1)
plt.xlabel("r (arcsec)")
plt.legend()
plt.xlim(0,120)
plt.show()
#fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/reports/report:27-08-2020/06-coneobs.png")
sys.exit()


