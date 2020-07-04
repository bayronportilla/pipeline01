import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
from photutils import EllipticalAnnulus
from mcmax3d_analysis.mcmax3d_observables import convert_flux
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_observation
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
import cflux_alma 
import cflux_jband
plt.style.use('fancy')

def output_data():
    
    ############################################################
    # Converting MCMax3D spectra 
    convert_flux('output/MCSpec0001.dat','spectrum_PDS70_system.dat')
    convert_flux('output/star0001.dat','spectrum_PDS70_star.dat')

    ############################################################
    # Loading surface density profile
    surface_density=np.loadtxt('surface_density_PDS70.dat')

    ############################################################
    # Loading data (converted) from MCMax3D 
    system_spectrum=np.loadtxt('spectrum_PDS70_system.dat')
    stellar_spectrum=np.loadtxt('spectrum_PDS70_star.dat')

    ############################################################
    # Loading photometric data
    photometric_data=np.loadtxt('PDS70_photometry.dat')
    x_photometric=np.reshape(photometric_data[:,0:1],photometric_data.shape[0])
    y_photometric=np.reshape(photometric_data[:,1:2],photometric_data.shape[0])
    y_photometric_error=np.reshape(photometric_data[:,2:3],photometric_data.shape[0])

    ############################################################
    # Creating arrays for plotting
    x_system=np.reshape(system_spectrum[:,0:1],system_spectrum.shape[0])
    y_system=np.reshape(system_spectrum[:,1:2],system_spectrum.shape[0])
    x_star=np.reshape(stellar_spectrum[:,0:1],stellar_spectrum.shape[0])
    y_star=np.reshape(stellar_spectrum[:,1:2],stellar_spectrum.shape[0])

    ############################################################
    # Creating array of errors
    y_system_error=x_photometric*y_photometric_error
   
    data_alma=cflux_alma.image("RTout0001_000854.89.fits.gz",0.074,0.057,63.0)
    data_jband=cflux_jband.image("RToutObs0001_000001.25.fits.gz")
    
    
    ############################################################
    # Dynamic range
    a_alma=0.01
    vmin_alma=np.percentile(data_alma,a_alma)
    vmax_alma=np.percentile(data_alma,100-a_alma)

    a_jband=0.35
    vmin_jband=np.percentile(data_jband,a_jband)
    vmax_jband=np.percentile(data_jband,100-a_jband)
    
  
    ############################################################
    # Observed profile

    odata=np.loadtxt("alma_radial_profile_observed.dat")
    r_obs=odata[:,0:1]
    b_obs=odata[:,1:2]

    odata_j_1=np.loadtxt("jband_radial_cut_0FWHM_smoothed.dat")
    r_obs_j_1=odata_j_1[:,0:1]
    b_obs_j_1=odata_j_1[:,1:2]


    ############################################################
    # Loading profiles modeled
    mprofile_alma=np.loadtxt("alma_radial_profile_modeled.dat")
    r_alma=np.reshape(mprofile_alma[:,0:1],mprofile_alma.shape[0])
    b_alma=np.reshape(mprofile_alma[:,1:2],mprofile_alma.shape[0])
    
    mprofile_jband=np.loadtxt("jband_radial_profile_modeled.dat")
    r_jband=np.reshape(mprofile_jband[:,0:1],mprofile_jband.shape[0])
    b_jband=np.reshape(mprofile_jband[:,1:2],mprofile_jband.shape[0])

    ############################################################
    # Plotting
    fig=plt.figure(figsize=(18,12))
    gs=gridspec.GridSpec(2,3)
    ax1=plt.subplot(gs[0,0:1])
    ax2=plt.subplot(gs[0,1:2])
    ax3=plt.subplot(gs[0,2:3])
    ax4=plt.subplot(gs[1,0:1])
    ax5=plt.subplot(gs[1,1:2])
    ax6=plt.subplot(gs[1,2:3])
    # Plot images
    path_fits_file='/output/RToutObs0001_000001.25.fits.gz'
    path_image_file='Image_jband.out'
    path_input_file='input.dat'
    # Loading Image.out info
    imfile=open(path_image_file).readlines()
    for line in imfile:
        if line.split('=')[0]=='MCobs:fov':
            fov=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:npix':
            npix=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:phi':
            phi=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:theta':
            theta=float(line.split('=')[1].split('!')[0])
        else:
            continue
    # Loading input.dat info
            

    # Determining extent
    #extent_angular=[0.5*fov,-0.5*fov,-0.5*fov,0.5*fov]
    xc=0#data_alma.shape[0]*0.5
    yc=0#data_alma.shape[0]*0.5
    
    infile=open(path_input_file).readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])
        elif line.split('=')[0]=='zone1:Rin':
            Rin01=float(line.split('=')[1])
        elif line.split('=')[0]=='zone1:Rout':
            Rout01=float(line.split('=')[1])
        elif line.split('=')[0]=='zone2:Rin':
            Rin02=float(line.split('=')[1])
        elif line.split('=')[0]=='zone2:Rout':
            Rout02=float(line.split('=')[1])
        elif line.split('=')[0]=='zone3:Rin':
            Rin03=float(line.split('=')[1])
        elif line.split('=')[0]=='zone3:Rout':
            Rout03=float(line.split('=')[1])
        else:
            continue
    
    a1=Rout01
    a2=Rout02
    

    e=np.sin((theta*units.deg).to(units.rad).value)
    d=(d*units.pc).value
    pxsize=fov/npix
 
    ax1.imshow(data_alma,clim=(vmin_alma,vmax_alma))
    ax2.imshow(data_jband,clim=(vmin_jband,vmax_jband))    
 
    msize=7.0
    lwidth=3.0
    ax3.plot(r_obs,b_obs,label="observation",linewidth=lwidth)
    ax3.plot(r_alma,b_alma,'.',markersize=msize,label="model")
    ax4.plot(r_obs_j_1,b_obs_j_1,"+",label="observation 0.0FWHM")
    ax4.plot(r_jband,b_jband,".",linewidth=lwidth,label="model")
    ax4.axvline(x=Rout01,linestyle='--',linewidth=0.5)
    ax4.axvline(x=Rout02,linestyle='--',linewidth=0.5)
#    ax4.axvline(x=Rout03,linestyle='--',linewidth=0.5)
    ax5.plot(x_system,x_system*y_system,label='MCMax3D output',linewidth=3.0)
    ax5.plot(x_star,x_star*y_star,color='orange')
    ax5.errorbar(x_photometric,x_photometric*y_photometric,yerr=y_system_error,label='Photometric data',fmt='.',markersize=10,color='green')
    ax6.plot(surface_density[:,0:1],surface_density[:,1:2])

    ############################################################
    # Scales
    ax5.set_xscale('log')
    ax5.set_yscale('log')
    #ax6.set_xscale('log')
    ax6.set_yscale('log')

    ############################################################
    # lims
    ax5.set_ylim(1e-17,1e-12)
    ax6.set_ylim(1e-5,1e1)
    ax4.set_ylim(-0.11,4.0)
    
    xmin,xmax=(200,800)
    ymin,ymax=(200,800)
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymax,ymin)
    ax2.set_xlim(xmin,xmax)
    ax2.set_ylim(ymax,ymin)
    

    ############################################################
    # Legend
    ax3.legend()
    ax4.legend(loc="upper left")
    
    ############################################################
    # Labels
    ax3.set_xlabel(r"$r$ (AU)")
    ax3.set_ylabel(r"Normalized density flux")
    ax4.set_xlabel(r"$r$ (AU)")
    ax4.set_ylabel(r"Normalized density flux")
    ax5.set_xlabel(r'$\lambda \, (\mu\mathrm{m})$')
    ax5.set_ylabel(r'$\lambda F_{\lambda}$ (W/m^2)')
    ax6.set_xlabel(r'$r$ (AU)')
    ax6.set_ylabel(r'$\Sigma_{\mathrm{dust}}$ (g/cm^2)')
    #plt.show()
    
    fig.savefig("fig_all.png")

    return None

#output_data()

