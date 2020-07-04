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
plt.style.use("fancy")

def shift(M):
    a=M.min()
    b=M.max()
    c=-1
    d=+1
    for i in range(0,M.shape[0]):
        for j in range(0,M.shape[1]):
            M[i][j]=c+((d-c)/(b-a))*(M[i][j]-a)
    return M


def make_plot(observation,pobs,pxsizeobs,model,pmod,pxsizemod,**kwargs):
    
    ############################################################
    #
    # Inputs.
    # observation: fits file of the observed image
    # pobs: peak flux value of the observed image at 54 AU
    # pxsizeobs: pixel scale of the observation in arcsex/px
    # model: fits file of the modeled image
    # pmod: peak flux value of the modeled image at 54 AU
    # pxsizemod: pixel scale of the model in arcsex/px
    #
    # Returns the final comparisson between the observed and 
    # modeled Qphi images. The compared images have the proper
    # orientation and the same color scale.
    #
    ############################################################


    ############################################################
    # Loading data
    hdu_obs=fits.open(observation)
    hdu_mod=fits.open(model)

    data_obs=hdu_obs[0].data/pobs
    data_mod=hdu_mod[0].data/pmod
    """
    plt.imshow(data_obs,clim=(-1,1))
    plt.show()
    plt.imshow(data_mod)
    plt.show()
    sys.exit()
    """

    ############################################################
    # Determine the fov
    fov_obs=data_obs.shape[0]*pxsizeobs
    fov_mod=data_mod.shape[0]*pxsizemod
    extent_obs=(0.5*fov_obs,-0.5*fov_obs,-0.5*fov_obs,0.5*fov_obs)
    extent_mod=(0.5*fov_mod,-0.5*fov_mod,-0.5*fov_mod,0.5*fov_mod)


    ############################################################
    # Limits of the final plot (those must be contained inside
    # extent_obs and extent_mod)
    lims=(+1.2,-1.2,-1.2,+1.2) # (xmax,xmin,ymin,ymax) arcsec
    
    
    ############################################################
    # General variables for plotting
    lsize=12 # Label size
    mapcolor="RdBu"


    ############################################################
    # Plotting
    fig,(ax_1,ax_2,cax)=plt.subplots(1,3,figsize=(8,4),
                                     gridspec_kw={"width_ratios":[1,1, 0.05]})                                     
    plt.subplots_adjust(wspace=0.0,hspace=0.0)

    # Limits on color scale
    a=0.1
    
    vmin_obs=-1
    vmax_obs=+1
    vmin_mod=-1
    vmax_mod=+1
    """
    vmin_obs=np.percentile(data_obs,a)
    vmax_obs=np.percentile(data_obs,100-a)
    vmin_mod=np.percentile(data_obs,a)
    vmax_mod=np.percentile(data_obs,100-a)
    """
    f_1=ax_1.imshow(data_obs,clim=(vmin_obs,vmax_obs),origin="lower",extent=extent_obs,cmap=mapcolor)
    f_2=ax_2.imshow(data_mod,clim=(vmin_mod,vmax_mod),origin="lower",extent=extent_mod,cmap=mapcolor)


    # Axes limits
    ax_1.set_xlim(lims[0],lims[1])
    ax_1.set_ylim(lims[2],lims[3])
    ax_2.set_xlim(lims[0],lims[1])
    ax_2.set_ylim(lims[2],lims[3])
    

    # Axes' ticks parameters
    ax_1.tick_params(top='on',right='on',labelright="off") 
    ax_2.tick_params(top='on',right='on',labelleft="off")
    ax_1.locator_params(axis='y',nbins=5)
    ax_2.locator_params(axis='y',nbins=5)

    
    # Color bar's parameters
    ip=InsetPosition(ax_2,[1,0,0.05,1]) # [(bar's left corner coordinates relative to ax_2),width,height]
    cax.set_axes_locator(ip)
    cbar=fig.colorbar(f_1,cax=cax)
    cbar.set_label(r"Normalized $Q_{\phi}$ signal",fontsize=lsize)

    
    # Common axes labels
    fig.text(0.5,0.0,r'$\Delta \mathrm{R.A.}$ (arcsec)', va='bottom', ha='center',fontsize=lsize)
    fig.text(0.06,0.5,r'$\Delta \mathrm{Dec.}$ (arcsec)', va='center', ha='center', rotation='vertical',fontsize=lsize)


    fig.savefig("output/Qphi.png")
    plt.show()
    
    return "File generated!"


B_max_obs=3.228573595978149
pxsize_obs=0.01226 

B_max_mod=0.136101396553597
pxsize_mod=0.004

make_plot("../../observations/PDS_70_2017-08-01_QPHI_amorph.fits",
          B_max_obs,pxsize_obs,
          "data/Qphi_model_rotated.fits",
          B_max_mod,pxsize_mod)



