import numpy as np
import matplotlib.pyplot as plt
from gofish import imagecube
import sys
import matplotlib.gridspec as gridspec


cube=imagecube('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',FOV=3)


"""
rvals, tvals, _ = cube.disk_coords(x0=0.0, y0=0.0, inc=49.7, PA=158.6)
fig, ax = plt.subplots()
im = ax.imshow(tvals, origin='lower', extent=cube.extent)
cb = plt.colorbar(im, pad=0.02)
ax.set_xlabel('Offset (arcsec)')
ax.set_ylabel('Offset (arcsec)')
cb.set_label('Polar Angle (radians)', rotation=270, labelpad=13)
plt.show()
sys.exit()
"""
"""
xm, ym, dym = cube.radial_profile(inc=49.7,PA=158.6,dist=113.43,x0=0.0,y0=0.0)
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
    #                               exclude_PA=False)

    return x,y,dy

padisk=158.6
incdisk=49.7
x0,y0,dy0=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                      padisk,incdisk,113.43,118.4,20,0.02)
x1,y1,dy1=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                      padisk,incdisk,113.43,90.0,20,0.02)
x2,y2,dy2=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                      padisk,incdisk,113.43,0.0,20,0.02)
x3,y3,dy3=get_profile('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits',
                      padisk,incdisk,113.43,-90.0,20,0.02)



x0=x0#*113.47
x1=x1#*113.47
x2=x2#*113.47
x3=x3#*113.47

fig=plt.figure(figsize=(5,12))
gs=gridspec.GridSpec(4,1,hspace=0)
ax0=plt.subplot(gs[0,0])
ax1=plt.subplot(gs[1,0])
ax2=plt.subplot(gs[2,0])
ax3=plt.subplot(gs[3,0])
ax0.errorbar(x0,(y0/max(ym)),yerr=dy0/max(ym),marker=".",fmt=" ",color="red",label="NW",capsize=1.25,capthick=1.25)
ax1.errorbar(x1,(y1/max(ym)),yerr=dy1/max(ym),marker=".",fmt=" ",color="green",label="SW",capsize=1.25,capthick=1.25)
ax2.errorbar(x2,(y2/max(ym)),yerr=dy2/max(ym),marker=".",fmt=" ",color="orange",label="SE",capsize=1.25,capthick=1.25)
ax3.errorbar(x3,(y3/max(ym)),yerr=dy3/max(ym),marker=".",fmt=" ",color="blue",label="NE",capsize=1.25,capthick=1.25)
ax0.legend()
ax1.legend()
ax2.legend()
ax3.legend()
ax0.axvline(0.65,0,2)
ax1.axvline(0.65,0,2)
ax2.axvline(0.65,0,2)
ax3.axvline(0.65,0,2)
maxlim=1.32
ax0.set_xlim(0,maxlim)
ax1.set_xlim(0,maxlim)
ax2.set_xlim(0,maxlim)
ax3.set_xlim(0,maxlim)


ax0.set_ylim(0,np.nanmax(y0/max(ym)))
ax1.set_ylim(0,np.nanmax(y0/max(ym)))
ax2.set_ylim(0,np.nanmax(y0/max(ym)))
ax3.set_ylim(0,np.nanmax(y0/max(ym)))


ax0.axes.get_xaxis().set_visible(False)
ax1.axes.get_xaxis().set_visible(False)
ax2.axes.get_xaxis().set_visible(False)

ax3.set_xlabel("$r$ (arcsec)")

#fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/reports/report:27-08-2020/shift.png")
plt.show()

sys.exit()
