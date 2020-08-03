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


def get_profile(file,pxsize,PA_aperture,PA_disk,inc,d,size,Nbins,**kwargs):

    ############################################################
    # 
    # file: the fits file of the observation
    # pxsize: pixel scale (arcsec/px)
    # PA_aperture: position angle of the aperture measured east-north (deg)
    # inc: disk's inclination (deg)
    # d: distance to the source (pc)
    # size: semi-major axis of the disk (AU)
    # 
    ############################################################


    ############################################################
    # Load data
    hdulist=fits.open(file)
    data_obs=hdulist[0].data[0][0]
    xc=hdulist[0].header['CRPIX1']
    yc=hdulist[0].header['CRPIX2']
    

    ############################################################
    # Derived properties
    angle_annulus=((PA_aperture-90.0)*units.deg).to(units.rad).value 
    e=np.sin((inc*units.deg).to(units.rad).value) 
    d_au=(d*units.pc).to(units.au).value # Distance (au)
    w=1.0
    h=1.0
    xc_array=[]
    yc_array=[]


    ############################################################
    # Creating elliptical aperture
    linear_lim=2*(size) # AU
    angular_lim=linear_lim/d_au # rad
    angular_lim=(angular_lim*units.rad).to(units.arcsec).value # arcsec
    pixel_lim=int(round(angular_lim/pxsize))
    dr=topx(30.0,pxsize,d) 
    a_in_array=[]

    for i in np.arange(yc+dr,yc+0.5*pixel_lim,dr):
        a_in_array.append(i-xc)
    a_out_array=[i+dr for i in a_in_array]
    b_out_array=[i*(1-e**2)**0.5 for i in a_out_array]

    apertures=[EllipticalAnnulus((yc,xc),a_in=ain,a_out=aout,b_out=bout,theta=angle_annulus)
               for (ain,aout,bout) in zip(a_in_array,a_out_array,b_out_array)]
        
    # Do a check?
    plt.imshow(data_obs)
    apertures[-1].plot(color='red',lw=1)
    plt.show()
        
    ############################################################
    # Define class "Bin"
    class Bin:
        def __init__(self,ID,theta_min,theta_max,plist):
            self.ID=ID
            self.theta_min=theta_min
            self.theta_max=theta_max
            self.plist=plist
        
        def getFlux(self):
            flux=0.0
            for pixel in self.plist:
                flux+=aperture_data[pixel[0],pixel[1]]
            return flux

        def getTheta(self):
            value=(self.theta_max-self.theta_min)*0.5+self.theta_min
            return value

    ############################################################
    # Creating array of bins
    bin_list=[]
    thetas=np.linspace(0,2*np.pi,Nbins+1)
    for i in range(0,Nbins):
        sbin=Bin(i+1,thetas[i],thetas[i+1],[])
        bin_list.append(sbin)

    
    M=np.zeros((Nbins,len(apertures)))

        
    for ii in range(0,len(apertures)):
        ############################################################
        # Creating aperture mask
        mask=apertures[ii].to_mask(method="center")
        """
        # Do a check?
        plt.imshow(mask)
        plt.colorbar()
        plt.show()
        """

        ############################################################
        # Extracting pixels located inside the aperture
        aperture_data=mask.multiply(data_obs)
        
        # Do a check?
        plt.imshow(aperture_data)
        plt.colorbar()
        plt.show()
           
        ############################################################
        # Creating array of pixel's index within the aperture 
        # relative to the star
        pixel_list=[]
        ycc=int(aperture_data.shape[0]*0.5)
        xcc=int(aperture_data.shape[1]*0.5)
        for i in range(0,aperture_data.shape[1]): # Over columns 
            for j in range(0,aperture_data.shape[0]): # Over rows
                if aperture_data[j,i]!=0.0:
                    pixel_list.append((j-ycc,i-xcc))


        ############################################################
        # Filling in bin_list
        for point in pixel_list:
            phi=np.arctan2(point[0],point[1])
            if phi<0.0:
                phi=2*np.pi+phi
            for sbin in bin_list:
                if sbin.theta_min<=phi<sbin.theta_max:
                    pixel=(int(point[0]+ycc),int(point[1]+xcc))
                    sbin.plist.append(pixel)


        ############################################################
        # Writing azimuthal profile 
        x=[]
        y=[]
        j=0
        for value in bin_list:
            M[j][ii]=value.getFlux()
            j+=1

    print(M)



get_profile("../PDS70/observations/PDS70_cont-final.fits",
            0.020,158.6,158.6,49.7,113.43,120.0,4,average=True)

