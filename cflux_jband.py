import numpy as np
import matplotlib.pyplot as plt
from astropy import units
from matplotlib.patches import Ellipse
from astropy.io import fits
from photutils import EllipticalAnnulus,CircularAnnulus,EllipticalAperture,RectangularAperture
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_observation
import sys
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

def combine_Polarizations(Q,U,phi0):
    Qphi= np.zeros(np.shape(Q))
    Uphi= np.zeros(np.shape(Q))
    x0= np.shape(Q)[0]*0.5-0.5
    y0= np.shape(Q)[1]*0.5-0.5
    phi0=(phi0*units.deg).to(units.rad).value
    for i in range(np.shape(Q)[0]): # over rows
            for j in range(np.shape(Q)[1]): # over columns
                phi= np.arctan2((float(i)-x0),(float(j)-y0))+phi0
                Qphi[i,j]= Q[i,j]*np.cos(2*phi)+U[i,j]*np.sin(2*phi)
                Uphi[i,j]= -Q[i,j]*np.sin(2*phi)+U[i,j]*np.cos(2*phi)
    return Qphi, Uphi

def shift(M,c,d):
    a=M.min()
    b=M.max()
    c=-1
    d=+1
    for i in range(0,M.shape[0]):
        for j in range(0,M.shape[1]):
            M[i][j]=c+((d-c)/(b-a))*(M[i][j]-a)
    return M


def image(fits_image):

    ############################################################
    # Absolute paths to files
    path_fits_image='../output/'+fits_image
    path_image_file='../Image_jband.out'
    path_input_file='../input.dat'

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
    phi=(phi_image*units.deg).to(units.rad).value # PA from north to east (rad)
    e=np.sin(theta) # eccentricity of the annulus

    ############################################################
    # Load MCMax3D image
    hdulist=fits.open(path_fits_image)
    data_Q=hdulist[0].data[1]
    data_U=hdulist[0].data[2]
    
    ############################################################
    # Start here
    Qphi_g, Uphi_g= combine_Polarizations(data_Q,data_U,0.0)
    
    # Shift data
    data_mod=Qphi_g

    return data_mod
    

def radial_profile(data,lim):

    ############################################################
    #
    # Get radial profiles - model
    #
    ############################################################

    ############################################################
    # Fetching information
    imfile=open("../Image_jband.out").readlines()
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

    infile=open("../input.dat").readlines()
    for line in infile:
        if line.split('=')[0]=='Distance':
            d=float(line.split('=')[1])

    ############################################################
    # Derived quantities
    pxsize=fov/npix # pixel scale (arcsec/px)
    phi=(phi_image*units.deg).to(units.rad).value # PA from north to east (rad)
    e=np.sin(theta) # eccentricity of the annulus

    angle_annulus=((0.0)*units.deg).to(units.rad).value 
    e=np.sin((theta*units.deg).to(units.rad).value) # eccentricity of the annulus

    # Determining limit for radial profile
    linear_lim=2*lim # AU
    angular_lim=linear_lim/((d*units.pc).to(units.au).value) # rad
    angular_lim=(angular_lim*units.rad).to(units.arcsec).value # arcsec
    pixel_lim=int(round(angular_lim/pxsize))

    xc=0.5*data.shape[0] # Image center in data coordinates
    yc=0.5*data.shape[1] # Image center in data coordinates
    dr=1.0 # Width of the annulus
    w=1.0
    h=1.0
    lim=120
    lim=toPX(lim,pxsize,d)
    InRad=np.pi/180.0
    xc_array=[]
    yc_array=[]

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
    InRad=np.pi/180.0
    apertures=RectangularAperture(positions,w,h,angle_annulus) 
    
    """
    ############################################################
    # Do a check
    a=0.4
    vmin=np.percentile(data,a)
    vmax=np.percentile(data,100-a)
    plt.imshow(data,clim=(vmin,vmax))
    plt.title("Image observation")
    apertures.plot(color='red',lw=1)
    #plt.xlim(400,600)
    #plt.ylim(600,400)
    plt.show()
    """
    phot_table=aperture_photometry(data,apertures)
    r_au=[toAU(phot_table['xcenter'][i].value,phot_table['ycenter'][i].value,xc,yc,pxsize,d) for i in range(0,len(phot_table))]
    brightness=[phot_table['aperture_sum'][i] for i in range(0,len(phot_table))]
    brightness=np.array(brightness)

    
    ############################################################
    # Creating brightness profile normalized
    rcmin=35.0
    rcmax=100.0
    bmaxc=[] 
    for i in range(0,len(r_au)):
        if rcmin<=r_au[i]<=rcmax:
            bmaxc.append(brightness[i])
    bmaxc=np.array(bmaxc)
    fac=1/max(bmaxc)
    brightness=brightness*fac
    
    """
    fig=plt.figure()
    ax=plt.axes()
    ax.plot(r_au,brightness,'.')
    ax.set_xlabel(r"Projected radial distance (AU)")
    ax.set_ylabel("Density flux (a.u.)")
    ax.set_title("Radial profile observation")
    #ax.set_ylim(-0.1,5)
    plt.show()
    sys.exit()
    """
    ############################################################
    # Creating file
    file=open('../jband_radial_profile_modeled.dat',"w")
    for i in range(0,len(r_au)):
        file.write('%.15e %.15e \n'%(r_au[i],brightness[i]))
    
    #data_mod=get_profile(data_mod,pxsize,lim)

    return None

#jband()
