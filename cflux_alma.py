import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_model
import sys
plt.style.use('fancy')

def image(fits_image,beam_x,beam_y,beam_angle):


    ############################################################
    # Absolute paths to files
    path_fits_image='output/'+fits_image
    path_image_file='Image_alma.out'
    path_input_file='input.dat'


    ############################################################
    # Fetching information
    imfile=open(path_image_file).readlines()
    for line in imfile:
        if line.split('=')[0]=='MCobs:fov':
            fov=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:npix':
            npix=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:phi':
            phi_image=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:theta':
            theta=float(line.split('=')[1].split('!')[0])
        else:
            continue

    infile=open(path_input_file).readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])


    ############################################################
    # Derived quantities
    pxsize=fov/npix # pixel scale (arcsec/px)
    theta=(theta*units.deg).to(units.rad).value # Inclination (rad)
    d=(d*units.pc).to(units.au).value # Distance (au)
    e=np.sin(theta) # eccentricity of the annulus


    ############################################################
    # Load MCMax3D image
    hdulist=fits.open(path_fits_image)
    data_mod=hdulist[0].data[0] # mJy/arcsec^2


    ############################################################
    # Convolve?
    #beam_x=0.074 # arcsec
    #beam_y=0.057 # arcsec
    data_mod=convolve_model(data_mod,fov,npix,beam_x,beam_y,beam_angle) # mJy/arcsec^2


    ############################################################
    # Convertion from mJy/arcsec^2 to mJy/beam
    beam_area=np.pi*(beam_x)*(beam_y)/(4*np.log(2)) 
    data_mod=data_mod*beam_area # mJy/beam
    print()
    print("Maximum value of the density flux in ALMA image (mJy/beam)",data_mod.max())
    print()

    return data_mod


def radial_profile(data,lim):

    ############################################################
    # Fetching information
    imfile=open("Image_alma.out").readlines()
    for line in imfile:
        if line.split('=')[0]=='MCobs:fov':
            fov=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:npix':
            npix=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:phi':
            phi_image=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:theta':
            theta=float(line.split('=')[1].split('!')[0])
        else:
            continue

    infile=open("input.dat").readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])


    ############################################################
    # Derived quantities
    pxsize=fov/npix # pixel scale (arcsec/px)
    theta=(theta*units.deg).to(units.rad).value # Inclination (rad)
    d=(d*units.pc).to(units.au).value # Distance (au)
    e=np.sin(theta) # eccentricity of the annulus


    ############################################################
    # Input params
    angle_annulus=0.0

    # Determining limit for radial profile
    #lim=120.0
    linear_lim=2*(lim) # AU
    angular_lim=linear_lim/d # rad
    angular_lim=(angular_lim*units.rad).to(units.arcsec).value # arcsec
    pixel_lim=int(round(angular_lim/pxsize))
    xc=0.5*data.shape[0] # Image center in data coordinates
    yc=0.5*data.shape[1] # Image center in data coordinates
    dr=1.0 # Width of the annulus
    a_in_array=[]
    for i in np.arange(yc+dr,yc+0.5*pixel_lim,dr):
        a_in_array.append(i-xc)
    a_out_array=[i+dr for i in a_in_array]
    b_out_array=[i*(1-e**2)**0.5 for i in a_out_array]

    apertures=[EllipticalAnnulus((yc,xc),a_in=ain,a_out=aout,b_out=bout,theta=angle_annulus)
               for (ain,aout,bout) in zip(a_in_array,a_out_array,b_out_array)]

    """
    ############################################################
    # Do a check
    a=0.01
    vmin_jband=np.percentile(data,a)
    vmax_jband=np.percentile(data,100-a)
    aperture=apertures[-1]
    plt.imshow(data,clim=(vmin_jband,vmax_jband))
    plt.title("Image model")
    aperture.plot(color='red',lw=1)
    plt.show()
    print(data.max())
    sys.exit()
    """

    # Radial distance of each annulus
    r_arcsec=[(j+0.5*(i-j))*pxsize for (i,j) in zip(a_out_array,a_in_array)] # arcsec
    r_rad=[(i*units.arcsec).to(units.rad).value for i in r_arcsec] # rad
    r_au=[(i*d) for i in r_rad] # AU

    # Creating numpy arrays
    r_au=np.array(r_au)
    r_arcsec=np.array(r_arcsec)

    phot_table=aperture_photometry(data,apertures)
    col_values=[]
    for col in phot_table.colnames:
        col_values.append(phot_table[col][0])
    brightness=[col_values[i] for i in range(3,len(col_values))]
    brightness=np.array(brightness)

    for i in range(0,len(brightness)):
        brightness[i]=brightness[i]/apertures[i].area

    rcmin=30.0
    rcmax=100.0
    bmaxc=[] 
    for i in range(0,len(r_au)):
        if rcmin<=r_au[i]<=rcmax:
            bmaxc.append(brightness[i])
    bmaxc=np.array(bmaxc)

    fac=1/max(bmaxc)
    brightness=brightness*fac

    """
    ############################################################
    # Creating brightness profile
    fig=plt.figure()
    ax=plt.axes()
    ax.plot(r_au,brightness/max(brightness),'*')
    ax.set_xlabel(r"Projected radial distance (AU)")
    ax.set_ylabel("Density flux (mJy/beam)")
    ax.set_title("Radial profile model")
    plt.show()
    sys.exit()
    """
    
    ############################################################
    # Creating file
    file=open('alma_radial_profile_modeled.dat',"w")
    for i in range(0,len(r_au)):
        file.write('%.15e %.15e \n'%(r_au[i],brightness[i]))

    return None

#alma()
