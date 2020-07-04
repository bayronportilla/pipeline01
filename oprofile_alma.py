import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture
from photutils import aperture_photometry
import sys
plt.style.use('fancy')


############################################################
# Loading observed image and reported parameters
file="../observations/PDS70_cont-final.fits"
pxsize=0.020 # arcsec/px
PA_disk=158.6 # deg
theta=49.7 # deg
d=113.43 # pc 
hdulist=fits.open(file)
data_obs=hdulist[0].data[0][0]*1000 # Converting to mJy/beam  


############################################################
# Input parameters
angle_annulus=((PA_disk-90.0)*units.deg).to(units.rad).value # Position angle
e=np.sin((theta*units.deg).to(units.rad).value) # eccentricity of the annulus


############################################################
# Determining boundaries for radial profile
linear_lim=2*(120.0) # AU
angular_lim=linear_lim/((d*units.pc).to(units.au).value) # rad
angular_lim=(angular_lim*units.rad).to(units.arcsec).value # arcsec
pixel_lim=int(round(angular_lim/pxsize)) # pixel


############################################################
# Set up aperture photometry
xc=0.5*data_obs.shape[0] # Image center in data coordinates
yc=0.5*data_obs.shape[1] # Image center in data coordinates
dr=1 # Width of the annulus
a_in_array=[]
for i in np.arange(yc+0.001,yc+0.5*pixel_lim+0.001,dr):
    a_in_array.append(i-xc)
a_out_array=[i+dr for i in a_in_array]
b_out_array=[i*(1-e**2)**0.5 for i in a_out_array]
apertures=[EllipticalAnnulus((yc,xc),a_in=ain,a_out=aout,b_out=bout,theta=angle_annulus)
           for (ain,aout,bout) in zip(a_in_array,a_out_array,b_out_array)]
"""
# Do a check
aperture=apertures[-1]
plt.imshow(data_obs)
aperture.plot(color='red',lw=1)
plt.show()
"""


# Radial distance of each annulus
r_arcsec=[(j+0.5*(i-j))*pxsize for (i,j) in zip(a_out_array,a_in_array)] # arcsec
r_rad=[(i*units.arcsec).to(units.rad).value for i in r_arcsec] # rad
r_au=[((i*d)*units.pc).to(units.au).value for i in r_rad] # AU


# Creating numpy arrays
r_au=np.array(r_au)
r_arcsec=np.array(r_arcsec)


############################################################
# Doing photometry
phot_table=aperture_photometry(data_obs,apertures)
col_values=[]
for col in phot_table.colnames:
    col_values.append(phot_table[col][0])
brightness=[col_values[i] for i in range(3,len(col_values))]
brightness=np.array(brightness)
for i in range(0,len(brightness)):
    brightness[i]=brightness[i]/apertures[i].area
"""
# Do a check?
fig=plt.figure()
ax=plt.axes()
ax.plot(r_au,brightness,'.')
ax.set_xlabel(r"Projected radial distance (AU)")
ax.set_ylabel("Density flux (mJy/beam)")
ax.set_title("Radial profile observation")
plt.show()
sys.exit()
"""


############################################################
# Finding position of the maximum value 
maxval=np.nanmax(data_obs)
jmax=np.where(data_obs==maxval)[0][0]
imax=np.where(data_obs==maxval)[1][0]
jc=xc
ic=yc
jmax=jmax-jc
imax=imax-ic
PA_max=(np.arctan2(jmax,imax)*units.rad).to(units.deg).value # w.r.t. x-axis
if PA_max<=0.0:
    PA_max=360.0+PA_max # w.r.t. x-axis and between 0 and 360 deg
PA_max=PA_max-90.0 # w.r.t. north 
if PA_max<=0.0:
    PA_max=360.0+PA_max # w.r.t north and between 0 and 360 deg
r_max=((((jmax**2+imax**2)**0.5)*pxsize)*units.arcsec).value # arcsec
r_max=d*r_max


############################################################
# Writing files
f1=open("info_max_alma.dat","w")
f2=open("alma_radial_profile_observed.dat","w")
f1.write("r_max=%.2f (AU)\n"%r_max)
f1.write("PA_max=%.2f (deg)\n"%PA_max)
f1.write("B_max=%.15f (mJy/beam)\n"%maxval)
for i in range(0,len(r_au)):
    f2.write("%.5e %.5e \n"%(r_au[i],brightness[i]/max(brightness)))
f1.close()
f2.close()

