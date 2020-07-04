from scipy import misc
import scipy
import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import sys
from mcmax3d_analysis.mcmax3d_observables import convert_flux,convert_flux_data
from mcmax3d_analysis.mcmax3d_image import display_image
from mcmax3d_analysis.mcmax3d_convolution import convolve_model
from astropy.convolution import Gaussian2DKernel


def pivot(data,r,PA,pxsize,d):

    ############################################################
    # Inputs parameters
    xc=data.shape[1]*0.5
    yc=data.shape[0]*0.5
    
    r=(r*units.au).to(units.pc).value # pc
    r=((r/d)*units.rad).to(units.arcsec).value # arcsec
    r=r/pxsize #px
    PA=((PA+90)*units.deg).to(units.rad).value # w.r.t. x-axis 

    xp=r*np.cos(PA)
    yp=r*np.sin(PA)

    x=xp+xc
    y=yp+yc

    return (x,y,data[int(round(y)),int(round(x))])
    


def prepare_Qphi_image(data,PA_disk):
    infile=open("input.dat").readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])

        
    ############################################################
    # Loading Image.out info
    imfile=open("Image_jband.out").readlines()
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

        
    ############################################################
    # Derived quantities
    pxsize=fov/npix # pixel scale (arcsec/px)
    phi=(phi*units.deg).to(units.rad).value # PA from north to east (rad)
    e=np.sin(theta) # eccentricity of the annulus

    
    ############################################################
    # Creating Qphi, Uphi
    data_rot=scipy.ndimage.rotate(data,-(PA_disk-90),reshape=False)


    ############################################################
    # Creating fits file
    hdu=fits.PrimaryHDU(data_rot)
    hdu.writeto("Qphi_model_rotated.fits",overwrite=True)

    return data_rot
        

def peak_flux_Qphi_model(Qphi_model_rotated):

    hdu=fits.open(Qphi_model_rotated)
    data_rot_Qphi=hdu[0].data

    infile=open("input.dat").readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])

        
    ############################################################
    # Loading Image.out info
    imfile=open("Image_jband.out").readlines()
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

        
    ############################################################
    # Derived quantities
    pxsize=fov/npix # pixel scale (arcsec/px)


    ############################################################
    # Searching for flux at a particular pixel 
    r_max_obs=54.64 # AU
    PA_max_obs=158.6 # deg
    xmax=pivot(data_rot_Qphi,r_max_obs,PA_max_obs,pxsize,d)[0]
    ymax=pivot(data_rot_Qphi,r_max_obs,PA_max_obs,pxsize,d)[1]
    Bmax=pivot(data_rot_Qphi,r_max_obs,PA_max_obs,pxsize,d)[2]

    
    return (xmax,ymax,pxsize,Bmax)


    
def prepare_alma_image(data,PA_disk,**kwargs):

    infile=open("input.dat").readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])

        
    ############################################################
    # Loading Image.out info
    imfile=open("Image_alma.out").readlines()
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

            
    ############################################################
    # Derived quantities
    pxsize=fov/npix # pixel scale (arcsec/px)
    phi=(phi*units.deg).to(units.rad).value # PA from north to east (rad)
    e=np.sin(theta) # eccentricity of the annulus


    ############################################################
    # Rotate image        
    data_rot=scipy.ndimage.rotate(data,-(PA_disk-90),reshape=False)


    # Creating fits file
    hdu=fits.PrimaryHDU(data_rot)
    hdu.writeto("alma_model_rotated.fits",overwrite=True)


    return data_rot


def peak_flux_alma_model(alma_model_rotated):

    hdu=fits.open(alma_model_rotated)
    data_rot_alma=hdu[0].data

    infile=open("input.dat").readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])

        
    ############################################################
    # Loading Image.out info
    imfile=open("Image_jband.out").readlines()
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

        
    ############################################################
    # Derived quantities
    pxsize=fov/npix # pixel scale (arcsec/px)


    ############################################################
    # Searching for flux at a particular pixel 
    r_max_obs=65.44 # AU
    PA_max_obs=326.31 # deg
    xmax=pivot(data_rot_alma,r_max_obs,PA_max_obs,pxsize,d)[0]
    ymax=pivot(data_rot_alma,r_max_obs,PA_max_obs,pxsize,d)[1]
    Bmax=pivot(data_rot_alma,r_max_obs,PA_max_obs,pxsize,d)[2]
        
    return (xmax,ymax,pxsize,Bmax)




