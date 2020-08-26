import numpy as np
import matplotlib.pyplot as plt
from gofish import imagecube
import sys
import matplotlib.gridspec as gridspec
from astropy.io import fits
from astropy.coordinates import Angle
import astropy.units as u

hdu_r=fits.open("/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits")
data=hdu_r[0].data[0][0]
hdr_r=hdu_r[0].header
print(hdr_r['ctype2'])

sys.exit()
"""
hdr_w=fits.Header()
hdr_w.append(('bunit',hdr_r['bunit'],None))
hdr_w.append(('bmaj',hdr_r['bmaj'],None))
hdr_w.append(('bmin',hdr_r['bmin'],None))
hdr_w.append(('bpa',hdr_r['bpa'],None))
hdr_w.append(('cdelt1',hdr_r['cdelt1'],None))
hdr_w.append(('cdelt2',hdr_r['cdelt2'],None))
hdr_w.append(('crpix1',hdr_r['crpix1'],None))
hdr_w.append(('crpix2',hdr_r['crpix2'],None))
hdr_w.append(('crval1',hdr_r['crval1'],None))
hdr_w.append(('crval2',hdr_r['crval2'],None))
#hdr_w.append(('restfreq',hdr_r['restfrq'],None))
hdr_w.append(('ctype1',hdr_r['ctype1'],None))
hdr_w.append(('ctype2',hdr_r['ctype2'],None))
hdu_w=fits.PrimaryHDU(data,header=hdr_w)
hdu_w.writeto("/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70.fits",overwrite=True)
sys.exit()
"""
#cube=imagecube('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',FOV=3)
cube2=imagecube('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70.fits',FOV=3)
#xm, ym, dym = cube.radial_profile(inc=49.7,PA=158.6,dist=113.43,x0=0.02,y0=0.02)
xm2, ym2, dym2 = cube2.radial_profile(inc=49.7,PA=158.6,dist=113.43,x0=0.02,y0=0.02)
#plt.errorbar(xm,ym,dym,fmt='.')
plt.errorbar(xm2,ym2,dym2,fmt='.')
plt.show()
sys.exit()
"""
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
                               assume_correlated=False,
                               x0=0.02,y0=0.02)
                               #exclude_PA=False)
    return x,y,dy

ddisk=113.43
drsample=0.02
padisk=158.6
incdisk=49.7
padir=118.4
#padir=-61.6
widir=20.0
nbins=3 # odd integer

fig=plt.figure()
for i in range(int(-0.5*(nbins-1)),int(0.5*(nbins-1))+1):
    angle=padir+i*widir
    x,y,dy=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                       padisk,incdisk,ddisk,angle,widir,drsample)
    plt.errorbar(x*ddisk,y/np.nanmax(ym),dy/np.nanmax(ym),fmt=".",markersize=4,elinewidth=1,
                 label="PA=%.1f deg"%(angle+padisk))
    

plt.axvspan(34.1,43.1,color='red',alpha=0.1)
plt.xlabel("r (arcsec)")
plt.legend()
plt.xlim(0,120)
plt.show()
#fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/reports/report:27-08-2020/06-coneobs.png")
sys.exit()


