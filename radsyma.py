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


def get_profile(file,pxsize,PA_disk,inc,d,size,Nbins,dr,**kwargs):

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

        xc=data_obs.shape[0]*0.5
        yc=data_obs.shape[1]*0.5
        """
        xc=hdulist[0].header['CRPIX1']
        yc=hdulist[0].header['CRPIX2']
        """
        

    elif kwargs['type']=='mod':
        hdulist=fits.open(file)
        data_obs=hdulist[0].data
        xc=data_obs.shape[0]*0.5
        yc=data_obs.shape[1]*0.5
    

    ############################################################
    # Derived properties
    angle_annulus=((PA_disk-90.0)*units.deg).to(units.rad).value 
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
            """
            if value<=0.5*np.pi:
                value=value+1.5*np.pi
            else:
                value=value-0.5*np.pi
            """
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
            

    #thetas=np.linspace(0,2*np.pi,Nbins+1)
    thetas=np.linspace(0,2*np.pi,Nbins+1)
    M=np.zeros((Nbins,len(apertures)))
    E_beam=np.zeros((Nbins,len(apertures)))
    E_pixel=np.zeros((Nbins,len(apertures)))
    

    a_in_array=[i*pxsize*d for i in a_in_array]
    a_out_array=[i*pxsize*d for i in a_out_array]
    a_mid=np.array([(j-i)*0.5+i for (j,i) in zip(a_out_array,a_in_array)])
    

    for ii in range(0,len(apertures)):

        ############################################################
        # Creating array of bins
        bin_list=[]

        for i in range(0,Nbins):
            sbin=Bin(i+1,thetas[i],thetas[i+1],[])
            bin_list.append(sbin)
        
        if ii==0:
            midtheta=[]
            midr=[]
            for value in bin_list:
                angle=value.getTheta()-0.5*np.pi
                if angle<0.0:
                    angle=2*np.pi+angle
                midtheta.append(angle)
            midtheta=np.array(midtheta)


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
        # Filling in bin_list
        for point in pixel_list:
            phi=np.arctan2(point[0],point[1])#-0.5*np.pi
            if phi<0.0:
                phi=2*np.pi+phi
            for sbin in bin_list:
                if sbin.theta_min<=phi<sbin.theta_max:
                    pixel=(int(point[0]+ycc),int(point[1]+xcc))
                    sbin.plist.append(pixel)


        ############################################################
        # Writing result
        j=0 # Count over position angles
        for value in bin_list:
            #value.showFlux()
            M[j][ii]=value.getFlux()/value.getArea(apertures[ii])
            E_beam[j][ii]=value.getError_beam(apertures[ii])
            E_pixel[j][ii]=value.getError_pixel(apertures[ii])
            j+=1
            #print()

    Mmax=M.max()
    print()
    print("Max value (Jy/beam/bin_area): %.1e"%(Mmax))
    print("Max error (per beam): %.1e"%(np.nanmax(E_beam)))
    print("Max error (per pixel): %.1e"%(np.nanmax(E_pixel)))

    for i in range(0,M.shape[0]):
        M[i]=M[i]/Mmax
        E_beam[i]=E_beam[i]/Mmax
        E_pixel[i]=E_pixel[i]/Mmax


    ############################################################
    # Plotting    
    index_to_swap=[]
    for i in range(0,len(midtheta)):
        if 1.5*np.pi<=midtheta[i]<2*np.pi:
            index_to_swap.append(i)
    
    M_new=np.zeros((M.shape[0],M.shape[1]))
    E_beam_new=np.zeros((E_beam.shape[0],E_beam.shape[1]))
    E_pixel_new=np.zeros((E_pixel.shape[0],E_pixel.shape[1]))
    
    M_new[M.shape[0]-len(index_to_swap):M.shape[0]]=M[0:len(index_to_swap)]
    M_new[0:M.shape[0]-len(index_to_swap)]=M[len(index_to_swap):M.shape[0]]
    M=M_new

    E_beam_new[E_beam.shape[0]-len(index_to_swap):E_beam.shape[0]]=E_beam[0:len(index_to_swap)]
    E_beam_new[0:E_beam.shape[0]-len(index_to_swap)]=E_beam[len(index_to_swap):E_beam.shape[0]]
    E_beam=E_beam_new

    E_pixel_new[E_pixel.shape[0]-len(index_to_swap):E_pixel.shape[0]]=E_pixel[0:len(index_to_swap)]
    E_pixel_new[0:E_pixel.shape[0]-len(index_to_swap)]=E_pixel[len(index_to_swap):E_pixel.shape[0]]
    E_pixel=E_pixel_new


    midtheta=np.sort(midtheta)


    fig=plt.figure(figsize=(5,12))
    gs=gridspec.GridSpec(int(Nbins*0.5),1,hspace=0)
    for i in range(0,int(Nbins*0.5)):
        ax=plt.subplot(gs[i,0])
        #       ax.plot(a_mid,np.reshape(M[i:i+1,:],M.shape[1]),'.',color="red")
        #       ax.plot(-a_mid,np.reshape(M[i+int(0.5*Nbins):i+1+int(0.5*Nbins),:],M.shape[1]),'.',color="red")
        ax.errorbar(a_mid,np.reshape(M[i:i+1,:],M.shape[1]),
                    yerr=np.reshape(E_beam[i:i+1,:],E_beam.shape[1]),marker=".",fmt="o",color="red")
        ax.errorbar(-a_mid,np.reshape(M[i+int(0.5*Nbins):i+1+int(0.5*Nbins),:],M.shape[1]),
                    yerr=np.reshape(E_beam[i+int(0.5*Nbins):i+1+int(0.5*Nbins),:],E_beam.shape[1]),
                    marker=".",fmt="o",color="red")
        """
        ax.errorbar(a_mid,np.reshape(M[i:i+1,:],M.shape[1]),
                    yerr=np.reshape(E_pixel[i:i+1,:],E_pixel.shape[1]),marker=".",fmt="o",color="red")
        ax.errorbar(-a_mid,np.reshape(M[i+int(0.5*Nbins):i+1+int(0.5*Nbins),:],M.shape[1]),
                    yerr=np.reshape(E_pixel[i+int(0.5*Nbins):i+1+int(0.5*Nbins),:],E_pixel.shape[1]),marker=".",fmt="o",color="red")
        """
        ax.axvline(+74,0,1)
        ax.axvline(-74,0,1)
        ax.tick_params(labelleft=False,left=False)
        ax.set_ylabel(r"%.1f"%((midtheta[i]*units.rad).to(units.deg).value))
        ax.set_xlabel(r"$r$(AU)")
    plt.show()


#get_profile("../PDS70/observations/PDS70_cont-final.fits",
#            0.020,158.6,49.7,113.43,120.0,18,4,type='obs')

get_profile("/data/users/bportilla/runs/final_runs/model_v.02.02/alma_model_rotated.fits",
            0.004,158.6,49.7,113.43,120.0,18,4,type='mod')
