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
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar
from mpl_toolkits.axes_grid1 import AxesGrid
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec
from matplotlib import patches
plt.style.use("fancy")



def make_plot(observation,maxobs,pxsizeobs,
              model,maxmod,pxsizemod,
              beamx,beamy,beam_angle,**kwargs):

    ############################################################
    #
    # Inputs.
    # observation: fits file of the observed image
    # pxsizeobs: pixel scale of the observation in arcsex/px
    # model: fits file of the modeled image
    # pxsizemod: pixel scale of the model in arcsex/px
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
    hdu_obs=fits.open(observation)
    hdu_mod=fits.open(model)

    data_obs=(hdu_obs[0].data[0][0])*1000.0 # mJy/beam
    data_mod=hdu_mod[0].data # mJy/beam

    data_obs=data_obs/maxobs
    data_mod=data_mod/maxmod


    ############################################################
    # Computing residuals
    data_res=data_obs-data_mod


    ############################################################
    # Determine the fov
    fov_obs=data_obs.shape[0]*pxsizeobs
    fov_mod=data_mod.shape[0]*pxsizemod
    extent_obs=(0.5*fov_obs,-0.5*fov_obs,-0.5*fov_obs,0.5*fov_obs)
    extent_mod=(0.5*fov_mod,-0.5*fov_mod,-0.5*fov_mod,0.5*fov_mod)


    ############################################################
    # Limits of the final plot (those must be contained inside
    # extent_obs and extent_mod)
    lims=(+1.25,-1.25,-1.25,+1.25) # (xmax,xmin,ymin,ymax) arcsec
    
    
    ############################################################
    # General variables for plotting
    lsize=12 # Label size
    mapcolor_images="inferno"
    mapcolor_residuals="coolwarm"
    
    
    ############################################################
    # Plotting
    
    fig=plt.figure(figsize=(8,3.6))

    gs = gridspec.GridSpec(2,3,hspace=0.0,wspace=0.0,height_ratios=[0.05,1])

    ax_1=plt.subplot(gs[0,0:2])
    ax_2=plt.subplot(gs[0,2])
    ax_3=plt.subplot(gs[1,0])
    ax_4=plt.subplot(gs[1,1])
    ax_5=plt.subplot(gs[1,2])

    ax_3.set_anchor("N")
    ax_4.set_anchor("N")
    ax_5.set_anchor("N")

    
    # Limits on color scale
    """
    a=0.01
    vmin_obs=np.percentile(data_obs,a)
    vmax_obs=np.percentile(data_obs,100-a)
    vmin_mod=np.percentile(data_obs,a)
    vmax_mod=np.percentile(data_obs,100-a)
    """
    if(abs(np.nanmin(data_res))>abs(np.nanmax(data_res))):
        absval=abs(np.nanmin(data_res))
    else:
        absval=abs(np.nanmax(data_res))
    
    f_1=ax_3.imshow(data_obs,origin="lower",extent=extent_obs,cmap=mapcolor_images)
    f_2=ax_4.imshow(data_mod,origin="lower",extent=extent_mod,cmap=mapcolor_images)
    f_3=ax_5.imshow(data_res,origin="lower",extent=extent_obs,cmap=mapcolor_residuals,clim=(-absval,absval))
    
    cbar_1=fig.colorbar(f_1,cax=ax_1,orientation="horizontal")
    cbar_2=fig.colorbar(f_3,cax=ax_2,orientation="horizontal")
    ax_1.xaxis.set_ticks_position("top")
    ax_2.xaxis.set_ticks_position("top")
    ax_1.xaxis.set_label_position("top")
    ax_2.xaxis.set_label_position("top")
    cbar_1.set_label(r"Normalized surface brightness")
    cbar_2.set_label(r"Residuals")


    # Axes limits
    ax_3.set_xlim(lims[0],lims[1])
    ax_3.set_ylim(lims[2],lims[3])
    ax_4.set_xlim(lims[0],lims[1])
    ax_4.set_ylim(lims[2],lims[3])
    ax_5.set_xlim(lims[0],lims[1])
    ax_5.set_ylim(lims[2],lims[3])
    
    
    # Axes' ticks parameters
    ax_3.tick_params(top='on',right='on',labelright="off") 
    ax_4.tick_params(top='on',right='on',labelleft="off",labelright="off")
    ax_5.tick_params(top='on',right='on',labelleft="off")
    #ax_3.spines['bottom'].set_color('red')
    #ax_3.tick_params(axis='x', colors='grey')
    ax_3.spines['right'].set_color('grey')
    ax_4.spines['left'].set_color('grey')
    ax_3.tick_params(which='major', color="grey")
    ax_4.tick_params(which='major', color="grey")
    
    ax_3.locator_params(axis='x',nbins=5)
    ax_4.locator_params(axis='x',nbins=5)
    ax_5.locator_params(axis='x',nbins=5)
    
    
    ############################################################
    # Draw beam
    e3=patches.Ellipse((1, -1),beamx,beamy, angle=beam_angle, linewidth=2, fill=True, zorder=2,color="darkgrey")
    e4=patches.Ellipse((1, -1),beamx,beamy, angle=beam_angle, linewidth=2, fill=True, zorder=2,color="darkgrey")
    e5=patches.Ellipse((1, -1),beamx,beamy, angle=beam_angle, linewidth=2, fill=True, zorder=2,color="black")
    ax_3.add_patch(e3)
    ax_4.add_patch(e4)
    ax_5.add_patch(e5)
    
    
    # Common axes labels
    fig.text(0.5,0.15,r'$\Delta \mathrm{R.A.}$ (arcsec)', va='bottom', ha='center',fontsize=lsize)
    fig.text(0.06,0.56,r'$\Delta \mathrm{Dec.}$ (arcsec)', va='center', ha='center', rotation='vertical',fontsize=lsize)
    
    #fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/reports/report:19-06-2020/%s.png"%(mapcolor))
    fig.savefig("output/residuals.png")
    plt.show()
    

    return "File generated!"

max_obs=1#1.702471971511841
pxsize_obs=0.020
beam_x=0.074
beam_y=0.057

max_mod=1#1.062021071191004
pxsize_mod=0.020



make_plot("../../observations/PDS70_cont-final.fits",max_obs,pxsize_obs,
          "data/alma_model_rotated.fits",max_mod,pxsize_mod,
          beam_x,beam_y,-153.0)



