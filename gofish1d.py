import numpy as np
import matplotlib.pyplot as plt
from gofish import imagecube
import sys
import matplotlib.gridspec as gridspec

"""
cube = imagecube('TWHya_CS_32.fits')
x, y, dy = cube.radial_profile(r_min=0.0, r_max=1.0,inc=5.0,
                               PA=152., mstar=0.88, dist=59.5,
                               PA_min=80,PA_max=100.0)

x2, y2, dy2 = cube.radial_profile(r_min=0.0, r_max=1.0, inc=5.0,
                               PA=152, mstar=0.88, dist=59.5,
                               PA_min=-100,PA_max=-80.0)

fig, ax = plt.subplots()
ax.errorbar(x, y, dy, fmt=' ', capsize=1.25, capthick=1.25, color='k', lw=1.0)
ax.errorbar(x2, y2, dy2, fmt=' ', capsize=1.25, capthick=1.25, color='k', lw=1.0)
ax.step(x, y, where='mid', color='k', lw=1.0)
ax.step(x2, y2, where='mid', color='k', lw=1.0)
plt.show()
sys.exit()
"""


def get_profile(file,PA_disk,inc,d,r_in,r_out,padir,widir,dr):
    cube=imagecube(file)
    x,y,dy=cube.radial_profile(inc=inc, PA=PA_disk, 
                               dist=d, 
                               PA_min=padir-0.5*widir,PA_max=padir+0.5*widir)

    return x,y,dy

x0,y0,dy0=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                      158.6,49.7,113.47,0.04,150.0,180.0,10,1)
x1,y1,dy1=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                      158.6,49.7,113.47,0.04,150.0,-90.0,10,1)
x2,y2,dy2=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                      158.6,49.7,113.47,0.04,150.0,0.0,10,1)
x3,y3,dy3=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                      158.6,49.7,113.47,0.04,150.0,90.0,10,1)



x0=x0*113.47
x1=x1*113.47
x2=x2*113.47
x3=x3*113.47
fig=plt.figure(figsize=(5,12))
gs=gridspec.GridSpec(4,1,hspace=0)
ax0=plt.subplot(gs[0,0])
ax1=plt.subplot(gs[1,0])
ax2=plt.subplot(gs[2,0])
ax3=plt.subplot(gs[3,0])
"""
ax0.errorbar(x0,y0/np.nanmax(y1),yerr=dy0/np.nanmax(y1),marker=".",fmt=" ",color="red",label="NW",capsize=1.25,elinewidth=0.5)
ax1.errorbar(x1,y1/np.nanmax(y1),yerr=dy1/np.nanmax(y1),marker=".",fmt=" ",color="green",label="SW",capsize=1.25,elinewidth=0.5)
ax2.errorbar(x2,y2/np.nanmax(y1),yerr=dy2/np.nanmax(y1),marker=".",fmt=" ",color="orange",label="SE",capsize=1.25,elinewidth=0.5)
ax3.errorbar(x3,y3/np.nanmax(y1),yerr=dy3/np.nanmax(y1),marker=".",fmt=" ",color="blue",label="NE",capsize=1.25,elinewidth=0.5)
"""
ax0.errorbar(x0,y0,yerr=dy0,marker=".",fmt=" ",color="red",label="NW",capsize=1.25,capthick=1.25)
ax1.errorbar(x1,y1,yerr=dy1,marker=".",fmt=" ",color="green",label="SW",capsize=1.25,capthick=1.25)
ax2.errorbar(x2,y2,yerr=dy2,marker=".",fmt=" ",color="orange",label="SE",capsize=1.25,capthick=1.25)
ax3.errorbar(x3,y3,yerr=dy3,marker=".",fmt=" ",color="blue",label="NE",capsize=1.25,capthick=1.25)
ax0.legend()
ax1.legend()
ax2.legend()
ax3.legend()

ax0.set_xlim(0,150)
ax1.set_xlim(0,150)
ax2.set_xlim(0,150)
ax3.set_xlim(0,150)

plt.show()

sys.exit()
