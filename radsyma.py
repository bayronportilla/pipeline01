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


def topx(l,pxsize,d):
    x=(l*units.au).to(units.pc).value
    x=((x/d)*units.rad).to(units.arcsec).value
    x=x/pxsize
    return x


def get_profile(file,pxsize,PA_aperture,PA_disk,inc,d,size,**kwargs):

    ############################################################
    # 
    # Get flux along the semi-major axis with a rectangular
    # aperture.
    #
    # file: the fits file of the observation
    # pxsize: pixel scale (arcsec/px)
    # PA_aperture: position angle of the aperture measured east-north (deg)
    # inc: disk's inclination (deg)
    # d: distance to the source (pc)
    # size: semi-major axis of the disk (AU)
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
    data_obs=hdulist[0].data[0][0]


    if kwargs['average']==False:
        ############################################################
        # Derived properties
        angle_annulus=((PA_aperture-90.0)*units.deg).to(units.rad).value 
        e=np.sin((inc*units.deg).to(units.rad).value) 
        xc=0.5*data_obs.shape[0] # Image center in data coordinates
        yc=0.5*data_obs.shape[1] # Image center in data coordinates
        w=1.0
        h=1.0
        lim=120.0
        lim=topx(lim,pxsize,d)
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
    
        # Do a check?
        a=0.01
        vmin=np.percentile(data_obs,a)
        vmax=np.percentile(data_obs,100-a)
        plt.imshow(data_obs)#,clim=(vmin,vmax))
        apertures.plot(color='red',lw=1)
        plt.show()
        #sys.exit()
    

        ############################################################
        # Performing aperture photometry
        phot_table=aperture_photometry(data_obs,apertures)
        r_au=[toAU(phot_table['xcenter'][i].value,phot_table['ycenter'][i].value,xc,yc,pxsize,d) for i in range(0,len(phot_table))]
        brightness=[phot_table['aperture_sum'][i] for i in range(0,len(phot_table))]
        for i in range(0,len(brightness)):
            brightness[i]=brightness[i]/apertures[i].area
    
        # Do a check?
        fig=plt.figure()
        ax=plt.axes()
        ax.plot(r_au,brightness,'.')
        ax.set_xlabel(r"Projected radial distance (AU)")
        ax.set_ylabel("Density flux (a.u.)")
        plt.show()
        #sys.exit()
    
        return r_au,brightness

    if kwargs['average']==True:

        ############################################################
        # Derived properties
        angle_annulus=((PA_aperture-90.0)*units.deg).to(units.rad).value 
        e=np.sin((inc*units.deg).to(units.rad).value) 
        xc=0.5*data_obs.shape[0] # Image center in data coordinates
        yc=0.5*data_obs.shape[1] # Image center in data coordinates
        w=1.0
        h=1.0
        lim=size
        lim=topx(lim,pxsize,d)
        xc_array=[]
        yc_array=[]


        ############################################################
        # Creating elliptical aperture
        a_disk=size
        b_disk=a_disk*(1-e**2)**0.5
        a_disk=topx(size,pxsize,d) 
        b_disk=topx(b_disk,pxsize,d)

        angle=((PA_disk+90)*units.deg).to(units.rad).value

        aperture=EllipticalAperture((xc,yc),a_disk,b_disk,theta=angle)
        
        
        # Do a check?
        plt.imshow(data_obs)
        aperture.plot(color='red',lw=1)
        plt.show()
        sys.exit()
        
    
        ############################################################
        # Creating aperture mask
        mask=aperture.to_mask(method="center")
        """
        # Do a check?
        plt.imshow(mask)
        plt.colorbar()
        plt.show()
        """
    
        ############################################################
        # Extracting pixels located inside the aperture
        aperture_data=mask.multiply(data)
        """
        # Do a check?
        plt.imshow(aperture_data)
        plt.colorbar()
        plt.show()
        """



get_profile("../PDS70/observations/PDS70_cont-final.fits",0.020,170.0,158.6,49.7,113.43,120.0,average=True)

