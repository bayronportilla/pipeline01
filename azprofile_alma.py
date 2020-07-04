import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_model
import sys
plt.style.use('fancy')

#Bmax_value=float(np.loadtxt("Bmax.dat"))

def topx(l,pxsize,d):
    x=(l*units.au).to(units.pc).value
    x=((x/d)*units.rad).to(units.arcsec).value
    x=x/pxsize
    return x


def azimuthal_profile(image,pxsize,amean,width,inc,PA,d,Nbins):

    ############################################################
    #
    # Inputs.
    # image: fits file of the image.
    # pxsize: pixel scale of the image in arcsex/px.
    # amean: mean semi-major axis of the elliptical aperture (AU).
    # width: total width of the aperture (AU).
    # inc: inclination of the disk (deg).
    # PA: east-north measured position angle of the disk (deg).
    # d: distance to the source in pc
    # Nbins: number of azimuthal bins 
    #
    # Returns 
    # The azimuthal brightness profile 
    # 
    #
    ############################################################

    
    ############################################################
    # Loading data
    """
    # Observation
    hdu=fits.open(image)
    data=hdu[0].data[0][0]*1000#/1.702471971511841
    """
    
    # Model
    hdu=fits.open(image)
    data=hdu[0].data#/Bmax_value
    
    
    ############################################################
    # Derived quantities
    x0=data.shape[0]*0.5
    y0=data.shape[1]*0.5
    inc=(inc*units.deg).to(units.rad).value
    e=np.sin(inc)


    ############################################################
    # Creating elliptical aperture
    amin=amean-width*0.5 # AU
    amax=amean+width*0.5 # AU
    bmax=amax*(1-e**2)**0.5 # AU

    amin=topx(amin,pxsize,d) # px
    amax=topx(amax,pxsize,d) # px
    bmax=topx(bmax,pxsize,d) # px

    angle=((PA+90)*units.deg).to(units.rad).value

    aperture=EllipticalAnnulus((x0,y0),a_in=amin,a_out=amax,b_out=bmax,theta=angle)
    
    """
    # Do a check?
    plt.imshow(data)
    aperture.plot(color='red',lw=1)
    plt.show()
    """
    
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

    
    ############################################################
    # Creating array of pixel's index within the aperture 
    # relative to the star
    pixel_list=[]
    yc=int(aperture_data.shape[0]*0.5)
    xc=int(aperture_data.shape[1]*0.5)
    for i in range(0,aperture_data.shape[1]): # Over columns 
        for j in range(0,aperture_data.shape[0]): # Over rows
            if aperture_data[j,i]!=0.0:
                pixel_list.append((j-yc,i-xc))
    

    ############################################################
    # Filling in bin_list
    for point in pixel_list:
        phi=np.arctan2(point[0],point[1])
        if phi<0.0:
            phi=2*np.pi+phi
        for sbin in bin_list:
            if sbin.theta_min<=phi<sbin.theta_max:
                pixel=(point[0]+yc,point[1]+xc)
                sbin.plist.append(pixel)
                

    ############################################################
    # Writing azimuthal profile 
    x=[]
    y=[]
    for value in bin_list:
        PA_bin=(value.getTheta()*units.rad).to(units.deg).value-90.0
        if PA_bin<0.0:
            PA_bin=360.0+PA_bin
        x.append(PA_bin)
        y.append(value.getFlux()/len(value.plist))
    
    f=open("../azprofile_alma_mod.dat","w")
    for i in range(0,len(x)):
        f.write("%.5f %.15f\n"%(x[i],y[i]))
    f.close()
    
    """
    plt.plot(x,y,".")
    plt.show()
    """
    
    return 0


