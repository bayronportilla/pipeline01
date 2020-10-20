from scipy import misc
import scipy
import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import patches
plt.style.use("fancy")

def pivot(data,pxsize):

    ############################################################
    # Inputs parameters
    xc=data.shape[1]*0.5
    yc=data.shape[0]*0.5
    d=113.43
    r=((54.64)*units.au).to(units.pc).value # pc
    r=((r/d)*units.rad).to(units.arcsec).value # arcsec
    r=r/pxsize #px
    PA=((158.6+90)*units.deg).to(units.rad).value # w.r.t. x-axis 
    
    
    xp=r*np.cos(PA)
    yp=r*np.sin(PA)

    x=xp+xc
    y=yp+yc

    return (x,y,data[int(round(y)),int(round(x))])
    

def shift(M):
    a=M.min()
    b=M.max()
    c=-1
    d=+1
    for i in range(0,M.shape[0]):
        for j in range(0,M.shape[1]):
            M[i][j]=c+((d-c)/(b-a))*(M[i][j]-a)
    return M


def make_plot(observation_alma,maxobs_alma,pxsizeobs_alma,model_alma,maxmod_alma,pxsizemod_alma,beamx,beamy,beam_angle,
              observation_Qphi,model_Qphi,
              **kwargs):

    ############################################################
    #
    # Inputs.
    # observation_alma: fits file of the observed image
    # pxsizeobs_alma: pixel scale of the observation in arcsex/px
    # model_alma: fits file of the modeled image
    # pxsizemod_alma: pixel scale of the model in arcsex/px
    # beamx: FWHM of the beam in the x direction (arcsec)
    # beamy: FWHM of the beam in the y direction (arcsec)
    # beam_angle: beam orientation respect to the x-axis
    #
    # Returns the final comparisson between the observed and 
    # modeled ALMA images. The residuals are also plotted. 
    #
    ############################################################


    ############################################################
    # Loading data
    hdu_obs_alma=fits.open(observation_alma)
    hdu_mod_alma=fits.open(model_alma)
    hdu_obs_Qphi=fits.open(observation_Qphi)
    hdu_mod_Qphi=fits.open(model_Qphi)

    data_obs_alma=(hdu_obs_alma[0].data[0][0])*1000.0 # mJy/beam
    data_mod_alma=hdu_mod_alma[0].data # mJy/beam
    data_obs_Qphi=hdu_obs_Qphi[0].data
    data_mod_Qphi=hdu_mod_Qphi[0].data

    data_obs_alma=data_obs_alma/maxobs_alma
    data_mod_alma=data_mod_alma/maxmod_alma

    ############################################################
    # Derive peak value of model
    imfile=open("../Image_jband.out").readlines()
    for line in imfile:
        if line.split('=')[0]=='MCobs:fov':
            fov_Qphi=float(line.split('=')[1].split('!')[0])
        elif line.split('=')[0]=='MCobs:npix':
            npix_Qphi=float(line.split('=')[1].split('!')[0])
        else:
            continue
    pxsizemod_Qphi=fov_Qphi/npix_Qphi # pixel scale model (arcsec/px)
    xmax_mod_Qphi=pivot(data_mod_Qphi,pxsizemod_Qphi)[0]
    ymax_mod_Qphi=pivot(data_mod_Qphi,pxsizemod_Qphi)[1]
    Bmax_mod_Qphi=pivot(data_mod_Qphi,pxsizemod_Qphi)[2]

    ############################################################
    # Constants of the observation
    Bmax_obs_Qphi=3.228573595978149
    pxsizeobs_Qphi=0.01226 


    data_obs_Qphi=data_obs_Qphi/Bmax_obs_Qphi
    data_mod_Qphi=data_mod_Qphi/Bmax_mod_Qphi

    ############################################################
    # Determine the fov
    fov_obs_alma=data_obs_alma.shape[0]*pxsizeobs_alma
    fov_mod_alma=data_mod_alma.shape[0]*pxsizemod_alma
    fov_obs_Qphi=data_obs_Qphi.shape[0]*pxsizeobs_Qphi
    fov_mod_Qphi=data_mod_Qphi.shape[0]*pxsizemod_Qphi

    extent_obs_alma=(0.5*fov_obs_alma,-0.5*fov_obs_alma,-0.5*fov_obs_alma,0.5*fov_obs_alma)
    extent_mod_alma=(0.5*fov_mod_alma,-0.5*fov_mod_alma,-0.5*fov_mod_alma,0.5*fov_mod_alma)
    extent_obs_Qphi=(0.5*fov_obs_Qphi,-0.5*fov_obs_Qphi,-0.5*fov_obs_Qphi,0.5*fov_obs_Qphi)
    extent_mod_Qphi=(0.5*fov_mod_Qphi,-0.5*fov_mod_Qphi,-0.5*fov_mod_Qphi,0.5*fov_mod_Qphi)

    ############################################################
    # Limits of the final plot (those must be contained inside
    # extent_obs and extent_mod)
    lims=(+1.25,-1.25,-1.25,+1.25) # (xmax,xmin,ymin,ymax) arcsec
    
    
    ############################################################
    # General variables for plotting
    lsize=14 # Label size
    mapcolor_alma="inferno"
    mapcolor_Qphi="RdBu"
    
    
    ############################################################
    # Plotting
    
    fig=plt.figure(figsize=(6.0,6.0))
    #gs = gridspec.GridSpec(2,3,hspace=0.0,wspace=0.0,width_ratios=[1,1,0.1],height_ratios=[1,1])
    gs = gridspec.GridSpec(2,3,hspace=0.0,wspace=0.0,width_ratios=[1.0,1.0,0.05])
    ax_1=plt.subplot(gs[0,0])
    ax_2=plt.subplot(gs[0,1])
    ax_3=plt.subplot(gs[0,2])
    ax_4=plt.subplot(gs[1,0])
    ax_5=plt.subplot(gs[1,1])
    ax_6=plt.subplot(gs[1,2])

    #ax_3.set_anchor("N")
    #ax_4.set_anchor("N")
    

    
    # Limits on color scale
    """
    a=0.01
    vmin_obs=np.percentile(data_obs_alma,a)
    vmax_obs=np.percentile(data_obs_alma,100-a)
    vmin_mod=np.percentile(data_obs_alma,a)
    vmax_mod=np.percentile(data_obs_alma,100-a)
    """
    # Limits on color scale
    

    f_1=ax_1.imshow(data_obs_alma,origin="lower",extent=extent_obs_alma,cmap=mapcolor_alma)
    f_2=ax_2.imshow(data_mod_alma,origin="lower",extent=extent_mod_alma,cmap=mapcolor_alma)
    f_4=ax_4.imshow(data_obs_Qphi,clim=(-1,1),origin="lower",extent=extent_obs_Qphi,cmap=mapcolor_Qphi)
    f_5=ax_5.imshow(data_mod_Qphi,clim=(-1,1),origin="lower",extent=extent_mod_Qphi,cmap=mapcolor_Qphi)
    
    
    cbar_1=fig.colorbar(f_1,cax=ax_3,orientation="vertical")
    cbar_2=fig.colorbar(f_5,cax=ax_6,orientation="vertical")
    #ax_3.xaxis.set_ticks_position("top")
    #ax_2.xaxis.set_ticks_position("top")
    #ax_3.xaxis.set_label_position("top")
    #ax_2.xaxis.set_label_position("top")
    cbar_1.set_label(r"Surface brightness (mJy/beam)",fontsize=10)
    cbar_2.set_label(r"Normalized $Q_{\phi}$ signal",fontsize=10)


    # Axes limits
    ax_1.set_xlim(lims[0],lims[1])
    ax_1.set_ylim(lims[2],lims[3])
    ax_2.set_xlim(lims[0],lims[1])
    ax_2.set_ylim(lims[2],lims[3])
    ax_4.set_xlim(lims[0],lims[1])
    ax_4.set_ylim(lims[2],lims[3])
    ax_5.set_xlim(lims[0],lims[1])
    ax_5.set_ylim(lims[2],lims[3])
    
    
    
    
    # Axes' ticks parameters
    ax_1.tick_params(top='on',right='on',labelright="off",labelbottom="off") 
    ax_2.tick_params(top='on',right='on',labelleft="off",labelright="off",labelbottom="off")
    ax_4.tick_params(top='on',right='on',labelright="off",labelsize=lsize) 
    ax_5.tick_params(top='on',right='on',labelleft="off",labelsize=lsize)
    ax_1.locator_params(axis='y',nbins=5)
    ax_2.locator_params(axis='y',nbins=5)
    ax_4.locator_params(axis='y',nbins=5)
    ax_5.locator_params(axis='y',nbins=5)
    ax_4.locator_params(axis='x',nbins=5)
    ax_5.locator_params(axis='x',nbins=5)
    
    #ax_3.spines['bottom'].set_color('red')
    #ax_3.tick_params(axis='x', colors='grey')
    ax_1.spines['right'].set_color('grey')
    ax_2.spines['left'].set_color('grey')
    ax_1.tick_params(which='major', color="grey")
    ax_2.tick_params(which='major', color="grey")
    
    ax_1.locator_params(axis='x',nbins=5)
    ax_2.locator_params(axis='x',nbins=5)
    
    
    ############################################################
    # Draw beam
    e3=patches.Ellipse((1, -1),beamx,beamy, angle=beam_angle, linewidth=2, fill=True, zorder=2,color="darkgrey")
    e4=patches.Ellipse((1, -1),beamx,beamy, angle=beam_angle, linewidth=2, fill=True, zorder=2,color="darkgrey")
    ax_1.add_patch(e3)
    ax_2.add_patch(e4)
    
        
    # Common axes labels
    fig.text(0.5,0.04,r'$\Delta \mathrm{R.A.}$ (arcsec)', va='bottom', ha='center',fontsize=lsize)
    fig.text(0.04,0.56,r'$\Delta \mathrm{Dec.}$ (arcsec)', va='center', ha='center', rotation='vertical',fontsize=lsize)
    
    name="Qphi_alma"
    fig.savefig("../%s.png"%(name))
    fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/paper_figures/%s.png"%(name),bbox_inches = 'tight')

    #fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/reports/report:19-06-2020/%s.png"%(mapcolor))
    #fig.savefig("output/residuals.png")
    plt.show()
    
    return "File generated!"

max_obs=1#1.702471971511841
pxsize_obs=0.020
beam_x=0.074
beam_y=0.057

max_mod=1#1.062021071191004
pxsize_mod=0.004



make_plot("/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS70_cont-final.fits",max_obs,pxsize_obs,
          "../alma_model_rotated.fits",max_mod,pxsize_mod,beam_x,beam_y,-153.0,
          "/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/PDS_70_2017-08-01_QPHI_amorph.fits",
          "../Qphi_model_rotated.fits")



