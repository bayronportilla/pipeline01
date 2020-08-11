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

    ############################################################
    #
    # d: distance to the source in pc
    # pxsize: pixel scale in arcsec/px
    #
    # Return: Length in AU
    #
    ############################################################

    dr=((x-xc)**2+(y-yc)**2)**0.5
    dr=(dr*pxsize)
    dr=(dr*units.arcsec).to(units.rad).value
    dr=dr*((d*units.pc).to(units.au).value)
    if x<=xc:
        return +1.0*dr
    else:
        return -1.0*dr


def topx(l,pxsize,d):

    ############################################################
    #
    # d: distance to the source in pc
    # pxsize: pixel scale in arcsec/px
    #
    # Return: Length in pixels
    #
    ############################################################

    x=(l*units.au).to(units.pc).value
    x=((x/d)*units.rad).to(units.arcsec).value
    x=x/pxsize
    return x


def get_profile(file,pxsize,PA_disk,inc,d,size,padir,widir,dr,**kwargs):
    
    ############################################################
    # 
    # Extract a radial cut of the brightness along diffent 
    # position angles.
    # 
    # file: the fits file of the observation
    # pxsize: pixel scale (arcsec/px)
    # PA_aperture: position angle of the aperture measured east-north (deg)
    # inc: disk's inclination (deg)
    # d: distance to the source (pc)
    # size: semi-major axis of the disk (AU)
    # padir: position angle of the desired direction (deg)
    # widir: width of the cone along the desired direction (deg)
    # dr: width of the annulus (AU)
    # 
    # The ouput is a matrix whose the rows and columns represent
    # the position angle and the radial distance of the flux 
    # measurement.
    #
    ############################################################


    ############################################################
    # Load ALMA data
    if kwargs['type']=='obs':
        hdulist=fits.open(file)
        data_obs=hdulist[0].data[0][0]
        xc=hdulist[0].header['CRPIX1']
        yc=hdulist[0].header['CRPIX2']

    elif kwargs['type']=='mod':
        hdulist=fits.open(file)
        data_obs=hdulist[0].data
        xc=data_obs.shape[0]*0.5
        yc=data_obs.shape[1]*0.5
    

    ############################################################
    # Derived properties
    angle_annulus=((PA_disk-90.0)*units.deg).to(units.rad).value 
    if padir<=270.0:
        padir=padir+90.0
    else:
        padir=padir-270.0
    e=np.sin((inc*units.deg).to(units.rad).value) 
    d_au=(d*units.pc).to(units.au).value 
    xc_array=[]
    yc_array=[]


    ############################################################
    # Creating elliptical aperture
    linear_lim=2*(size) # AU
    angular_lim=linear_lim/d_au # rad
    angular_lim=(angular_lim*units.rad).to(units.arcsec).value # arcsec
    pixel_lim=int(round(angular_lim/pxsize))
    dr=topx(dr,pxsize,d) # width of each annular aperture
    a_in_array=[]
    for i in np.arange(yc+dr,yc+0.5*pixel_lim,dr):
        a_in_array.append(i-xc)
    a_out_array=[i+dr for i in a_in_array]
    b_out_array=[i*(1-e**2)**0.5 for i in a_out_array]
    a_in_array=np.array(a_in_array)
    a_out_array=np.array(a_out_array)
    apertures=[EllipticalAnnulus((yc,xc),a_in=ain,a_out=aout,b_out=bout,theta=angle_annulus)
               for (ain,aout,bout) in zip(a_in_array,a_out_array,b_out_array)]

    print("Number of annular apertures: %d"%len(apertures))
    """
    # Do a check?
    plt.imshow(data_obs)
    apertures[-1].plot(color='red',lw=1)
    plt.show()
    """

    ############################################################
    # Determine Nbins
    Nbins=round(360.0/widir)


    ############################################################
    # Define class "Bin"
    class Bin:
        def __init__(self,ID,theta_min,theta_max,plist):
            self.ID=ID
            self.theta_min=theta_min
            self.theta_max=theta_max
            self.plist=plist

        def showFlux(self):
            i=0
            for pixel in self.plist:
                print(i,aperture_data[pixel[0],pixel[1]])
                i+=1

        def getArea(self,aperture):
            value=aperture.area/Nbins
            return value

        def getFlux(self):
            flux=0.0
            for pixel in self.plist:
                flux+=aperture_data[pixel[0],pixel[1]]
            return flux

        def getTheta(self):
            value=(self.theta_max-self.theta_min)*0.5+self.theta_min
            return value

        def getError_beam(self,aperture):
            beam_x=0.074 # arcsec
            beam_y=0.057 # arcsec
            flux_array=[]
            for pixel in self.plist:
                flux_array.append(aperture_data[pixel[0],pixel[1]])
            area=aperture.area/Nbins
            flux_array=np.array(flux_array)/area
            sigma=np.std(flux_array)
            beam_area=np.pi*(beam_x)*(beam_y)/(4*np.log(2)) 
            Nbeam=((aperture.area*pxsize**2)/Nbins)/beam_area
            return sigma/(Nbeam)**0.5

        def getError_pixel(self,aperture):
            flux_array=[]
            for pixel in self.plist:
                flux_array.append(aperture_data[pixel[0],pixel[1]])
            area=aperture.area/Nbins
            flux_array=np.array(flux_array)/area
            sigma=np.std(flux_array)
            Npixel=len(self.plist)
            return sigma/(Npixel)**0.5
            
    M=[]
    E_beam=[]
    E_pixel=[]
    
    a_in_array=[i*pxsize*d for i in a_in_array]
    a_out_array=[i*pxsize*d for i in a_out_array]
    a_mid=np.array([(j-i)*0.5+i for (j,i) in zip(a_out_array,a_in_array)])
    

    for ii in range(0,len(apertures)):

        ############################################################
        # Creating bin
        sbin=Bin(ii,padir-0.5*widir,padir+0.5*widir,[])        


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
        """
        # Do a check?
        plt.imshow(aperture_data)
        plt.colorbar()
        plt.show()
        """

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
        # Filling in sbin data
        for point in pixel_list:
            phi=np.arctan2(point[0],point[1])
            if phi<0.0:
                phi=2*np.pi+phi
            phi=phi*180.0/np.pi
            if sbin.theta_min<=phi<sbin.theta_max:
                pixel=(int(point[0]+ycc),int(point[1]+xcc))
                sbin.plist.append(pixel)


        ############################################################
        # Writing result
        #value.showFlux()
        M.append(sbin.getFlux()/sbin.getArea(apertures[ii]))
        E_beam.append(sbin.getError_beam(apertures[ii]))
        E_pixel.append(sbin.getError_pixel(apertures[ii]))

    M=np.array(M)
    E_beam=np.array(E_beam)
    E_pixel=np.array(E_pixel)
    """
    print()
    print("Max value (Jy/beam/bin_area): %.1e"%(max(M)))
    print("Max error (per beam): %.1e"%(np.nanmax(E_beam)))
    print("Max error (per pixel): %.1e"%(np.nanmax(E_pixel)))
    """
    
    """
    ############################################################
    # Plotting
    fig=plt.figure(figsize=(5,12))
    ax=plt.axes()
    ax.errorbar(a_mid,M,yerr=E_beam,marker=".",fmt="-",color="red")
    #   ax.tick_params(labelleft=False,left=False)
    #    ax.set_ylabel(r"%.1f"%((midtheta[i]*units.rad).to(units.deg).value))
    ax.set_xlabel(r"$r$(AU)")
    plt.show()
    """
    return a_mid,M,E_beam,E_pixel


"""
x1=get_profile("../PDS70/observations/PDS70_cont-final.fits",
               0.020,158.6,49.7,113.43,120.0,277,20,4,type='obs')[0]
y1=get_profile("../PDS70/observations/PDS70_cont-final.fits",
               0.020,158.6,49.7,113.43,120.0,257,20,4,type='obs')[1]
y2=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,277,20,4,type='obs')[1]
y3=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,297,20,4,type='obs')[1]
y4=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,237,20,4,type='obs')[1]
y5=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,317,20,4,type='obs')[1]
e1=get_profile("../PDS70/observations/PDS70_cont-final.fits",
               0.020,158.6,49.7,113.43,120.0,257,20,4,type='obs')[2]
e2=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,277,20,4,type='obs')[2]
e3=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,297,20,4,type='obs')[2]
e4=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,237,20,4,type='obs')[2]
e5=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,317,20,4,type='obs')[2]

plt.errorbar(x1,y4/max(y4),yerr=e4/max(y4),marker=".",fmt="--",color="orange",label="$237^\circ$")
plt.errorbar(x1,y1/max(y1),yerr=e1/max(y1),marker=".",fmt="--",color="red",label="$257^\circ$")
plt.errorbar(x1,y2/max(y2),yerr=e2/max(y2),marker=".",fmt="-",color="blue",label="$277^\circ$")
plt.errorbar(x1,y3/max(y3),yerr=e3/max(y3),marker=".",fmt="--",color="green",label="$297^\circ$")
plt.errorbar(x1,y5/max(y5),yerr=e5/max(y5),marker=".",fmt="--",color="magenta",label=r"$317^\circ$")
plt.xlabel(r"$r$ (AU)")
plt.axvspan(34.1,43.1, alpha=0.2, color='red')
plt.legend()
plt.show()
"""

x1=get_profile("../PDS70/observations/PDS70_cont-final.fits",
               0.020,158.6,49.7,113.43,120.0,57,20,4,type='obs')[0]
y1=get_profile("../PDS70/observations/PDS70_cont-final.fits",
               0.020,158.6,49.7,113.43,120.0,57,20,4,type='obs')[1]
y2=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,77,20,4,type='obs')[1]
y3=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,97,20,4,type='obs')[1]
y4=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,117,20,4,type='obs')[1]
y5=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,137,20,4,type='obs')[1]

e1=get_profile("../PDS70/observations/PDS70_cont-final.fits",
               0.020,158.6,49.7,113.43,120.0,57,20,4,type='obs')[2]
e2=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,77,20,4,type='obs')[2]
e3=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,97,20,4,type='obs')[2]
e4=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,117,20,4,type='obs')[2]
e5=get_profile("../PDS70/observations/PDS70_cont-final.fits",
                  0.020,158.6,49.7,113.43,120.0,137,20,4,type='obs')[2]

plt.errorbar(x1,y1/max(y1),yerr=e1/max(y1),marker=".",fmt="--",color="red",label="$57^\circ$")
plt.errorbar(x1,y2/max(y2),yerr=e2/max(y2),marker=".",fmt="-",color="blue",label="$77^\circ$")
plt.errorbar(x1,y3/max(y3),yerr=e3/max(y3),marker=".",fmt="--",color="green",label="$97^\circ$")
plt.errorbar(x1,y4/max(y4),yerr=e4/max(y4),marker=".",fmt="--",color="orange",label="$117^\circ$")
plt.errorbar(x1,y5/max(y5),yerr=e5/max(y5),marker=".",fmt="--",color="magenta",label=r"$137^\circ$")
plt.xlabel(r"$r$ (AU)")
plt.legend()
plt.show()


sys.exit()
get_profile("../PDS70/observations/PDS70_cont-final.fits",
               0.020,158.6,49.7,113.43,120.0,277,20,4,type='obs')
get_profile("/data/users/bportilla/runs/final_runs/model_v.03/alma_model_rotated.fits",
               0.004,158.6,49.7,113.43,120.0,297,20,4,type='mod')
