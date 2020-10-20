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
# Inputs
hdu_mod_1=fits.open("/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run02/alma_model_rotated.fits")
hdu_mod_2=fits.open("/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run03/alma_model_rotated.fits")
hdu_mod_3=fits.open("/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run04/alma_model_rotated.fits")
pxsizemod_1=0.004
pxsizemod_2=0.004
pxsizemod_3=0.004
beamx=0.074
beamy=0.057
beam_angle=-153.0

############################################################
# Loading data
data_mod_1=hdu_mod_1[0].data # mJy/beam
data_mod_2=hdu_mod_2[0].data # mJy/beam
data_mod_3=hdu_mod_3[0].data # mJy/beam


############################################################
# Determine pxsize
fov_mod_1=data_mod_1.shape[0]*pxsizemod_1
fov_mod_2=data_mod_2.shape[0]*pxsizemod_2
fov_mod_3=data_mod_3.shape[0]*pxsizemod_3
extent_mod_1=(0.5*fov_mod_1,-0.5*fov_mod_1,-0.5*fov_mod_1,0.5*fov_mod_1)
extent_mod_2=(0.5*fov_mod_2,-0.5*fov_mod_2,-0.5*fov_mod_2,0.5*fov_mod_2)
extent_mod_3=(0.5*fov_mod_3,-0.5*fov_mod_3,-0.5*fov_mod_3,0.5*fov_mod_3)


############################################################
# Limits of the final plot (those must be contained inside
# extent_obs and extent_mod)
lims=(+1.25,-1.25,-1.25,+1.25) # (xmax,xmin,ymin,ymax) arcsec

    
############################################################
# General variables for plotting
lsize=12 # Label size
mapcolor_images="inferno"

        
############################################################
# Plotting

fig=plt.figure(figsize=(8,3.6))
#fig=plt.figure(figsize=(5,4))
gs = gridspec.GridSpec(2,3,hspace=0.0,wspace=0.0,height_ratios=[0.05,1])

ax_1=plt.subplot(gs[0,0:3])
ax_2=plt.subplot(gs[1,0])
ax_3=plt.subplot(gs[1,1])
ax_4=plt.subplot(gs[1,2])

ax_1.set_anchor("S")
ax_2.set_anchor("N")
ax_3.set_anchor("N")
ax_4.set_anchor("N")


# Limits on color scale
"""
a=0.01
vmin_obs=np.percentile(data_obs,a)
vmax_obs=np.percentile(data_obs,100-a)
vmin_mod=np.percentile(data_obs,a)
vmax_mod=np.percentile(data_obs,100-a)
"""

f_2=ax_2.imshow(data_mod_1,origin="lower",extent=extent_mod_1,cmap=mapcolor_images)
f_3=ax_3.imshow(data_mod_2,origin="lower",extent=extent_mod_2,cmap=mapcolor_images)
f_4=ax_4.imshow(data_mod_3,origin="lower",extent=extent_mod_3,cmap=mapcolor_images)

cbar_1=fig.colorbar(f_2,cax=ax_1,orientation="horizontal")
cbar_1.ax.xaxis.set_tick_params(color='white')
ax_1.xaxis.set_ticks_position("top")
ax_1.xaxis.set_label_position("top")
cbar_1.set_label(r"Surface brightness (mJy/beam)")


# Axes limits
ax_2.set_xlim(lims[0],lims[1])
ax_2.set_ylim(lims[2],lims[3])
ax_3.set_xlim(lims[0],lims[1])
ax_3.set_ylim(lims[2],lims[3])
ax_4.set_xlim(lims[0],lims[1])
ax_4.set_ylim(lims[2],lims[3])
    
   
# Axes' ticks parameters
ax_2.tick_params(top='on',right='on',labelright="off") 
ax_3.tick_params(top='on',right='on',labelleft="off",labelright="off")
ax_4.tick_params(top='on',right='on',labelleft="off")

cspines="grey"
ax_2.spines['left'].set_color(cspines)
ax_2.spines['right'].set_color(cspines)
ax_2.spines['top'].set_color(cspines)
ax_2.spines['bottom'].set_color(cspines)

ax_3.spines['left'].set_color(cspines)
ax_3.spines['right'].set_color(cspines)
ax_3.spines['top'].set_color(cspines)
ax_3.spines['bottom'].set_color(cspines)

ax_4.spines['left'].set_color(cspines)
ax_4.spines['right'].set_color(cspines)
ax_4.spines['top'].set_color(cspines)
ax_4.spines['bottom'].set_color(cspines)

ax_2.tick_params(which='major', color=cspines)
ax_3.tick_params(which='major', color=cspines)
ax_4.tick_params(which='major', color=cspines)

ax_2.locator_params(axis='x',nbins=5)
ax_3.locator_params(axis='x',nbins=5)
ax_4.locator_params(axis='x',nbins=5)


############################################################
# Draw beam
e2=patches.Ellipse((1, -1),beamx,beamy, angle=beam_angle, linewidth=2, fill=True, zorder=2,color="darkgrey")
e3=patches.Ellipse((1, -1),beamx,beamy, angle=beam_angle, linewidth=2, fill=True, zorder=2,color="darkgrey")
e4=patches.Ellipse((1, -1),beamx,beamy, angle=beam_angle, linewidth=2, fill=True, zorder=2,color="darkgrey")
ax_2.add_patch(e2)
ax_3.add_patch(e3)
ax_4.add_patch(e4)

ax_3.set_xlabel(r'$\Delta \mathrm{R.A.}$ (arcsec)') 
ax_2.set_ylabel(r'$\Delta \mathrm{Dec.}$ (arcsec)')

############################################################
# Add planet mass info
tcolor="white"
ax_2.text(1,1,r"$0.5 \, M_\mathrm{Jup}$",color=tcolor)
ax_3.text(1,1,r"$5 \, M_\mathrm{Jup}$",color=tcolor)
ax_4.text(1,1,r"$12 \, M_\mathrm{Jup}$",color=tcolor)




plt.tight_layout()

fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/paper_figures/CPD_submm.png")


plt.show()
sys.exit()



