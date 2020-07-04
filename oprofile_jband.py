import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture,RectangularAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_observation
import matplotlib.gridspec as gridspec
import sys
from astropy.table import Table
plt.style.use('fancy')


def toAU(x,y,xc,yc,pxsize,d):
    dr=((x-xc)**2+(y-yc)**2)**0.5
    dr=(dr*pxsize)
    dr=(dr*units.arcsec).to(units.rad).value
    dr=dr*((d*units.pc).to(units.au).value)
    if x<=xc:
        return +1.0*dr
    else:
        return -1.0*dr


def toPX(l,pxsize,d):
    dl=(l*units.au).to(units.pc).value
    dl=((dl/d)*units.rad).to(units.arcsec).value
    dl=dl/pxsize
    return dl


def get_profile(file,pxsize,PA_disk,inc,d):

    ############################################################
    # 
    # Get flux along the semi-major axis with a rectangular
    # aperture.
    #
    # file: the fits file of the observation
    # pxsize: pixel scale (arcsec/px)
    # PA_disk: position angle of the disk measured east-north (deg)
    # inc: disk's inclination (deg)
    # d: distance to the source (pc)
    # 
    # returns a file with the brightness profile along the 
    # semi-major axis and a file with the position and flux 
    # value of the pixel located along the semi-major axis
    # in the "South-East" quadrant containing the maximum flux
    # value. The brightness profile returned is not normalized.
    # 
    ############################################################


    ############################################################
    # Load data
    hdulist=fits.open(file)
    data_obs=hdulist[0].data # adu's


    ############################################################
    # Derived properties
    angle_annulus=((PA_disk-90.0)*units.deg).to(units.rad).value 
    e=np.sin((inc*units.deg).to(units.rad).value) 
    xc=0.5*data_obs.shape[0] # Image center in data coordinates
    yc=0.5*data_obs.shape[1] # Image center in data coordinates
    w=1.0
    h=1.0
    lim=120.0
    lim=toPX(lim,pxsize,d)
    xc_array=[]
    yc_array=[]


    ############################################################
    # Setting up aperture photometry
    x0=xc
    y0=yc
    xc_array.append(x0)
    yc_array.append(y0)
    for l in np.arange(-lim,lim,1.0):
        xval=x0+l*np.cos(angle_annulus)
        yval=y0+l*np.sin(angle_annulus)
        xc_array.append(xval)
        yc_array.append(yval)
    positions=[(i,j) for (i,j) in zip(xc_array,yc_array)]
    apertures=RectangularAperture(positions,w,h,angle_annulus) 
    """
    # Do a check?
    a=0.01
    vmin=np.percentile(data_obs,a)
    vmax=np.percentile(data_obs,100-a)
    plt.imshow(data_obs,clim=(vmin,vmax))
    apertures.plot(color='red',lw=1)
    plt.show()
    sys.exit()
    """


    ############################################################
    # Performing aperture photometry
    phot_table=aperture_photometry(data_obs,apertures)
    r_au=[toAU(phot_table['xcenter'][i].value,phot_table['ycenter'][i].value,xc,yc,pxsize,d) for i in range(0,len(phot_table))]
    brightness=[phot_table['aperture_sum'][i] for i in range(0,len(phot_table))]
    for i in range(0,len(brightness)):
        brightness[i]=brightness[i]/apertures[i].area
    """
    # Do a check?
    fig=plt.figure()
    ax=plt.axes()
    ax.plot(r_au,brightness,'.')
    ax.set_xlabel(r"Projected radial distance (AU)")
    ax.set_ylabel("Density flux (a.u.)")
    plt.show()
    sys.exit()
    """


    ############################################################
    # Finding maximum value along semi-major axis
    rstart=40.0
    rend=80.0
    bmax=[] 
    for i in range(0,len(r_au)):
        if rstart<=r_au[i]<=rend:
            bmax.append(brightness[i])
    bmax=np.array(bmax)
    imax=np.where(brightness==max(bmax))[0][0]


    ############################################################
    # Writing files
    f1=open("info_max_Qphi.dat","w")
    f2=open("radial_cut_0FWHM.dat","w")
    f1.write("r_max=%.2f (AU)\n"%r_au[imax])
    f1.write("PA_max=%.2f (deg)\n"%PA_disk)
    f1.write("B_max=%.15f (a.u.)\n"%brightness[imax])
    for i in range(0,len(r_au)):
        f2.write("%.5e %.5e \n"%(r_au[i],brightness[i]))
    f1.close()
    f2.close()
    

    return None


get_profile("../observations/PDS_70_2017-08-01_QPHI_amorph.fits",0.01226,158.6,49.7,113.43)
