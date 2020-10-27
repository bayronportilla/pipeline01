import numpy as np
import sys
from astropy import units
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys
from photutils import EllipticalAnnulus
from mcmax3d_analysis.mcmax3d_observables import convert_flux
from photutils import aperture_photometry
from mcmax3d_analysis.mcmax3d_convolution import convolve_observation
plt.style.use('fancy')


#sim='/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run01'
sim='/data/users/bportilla/runs/final_runs/pds70'


############################################################
# Loading surface density profile
surface_density=np.loadtxt(sim+'/surface_density_PDS70.dat')

############################################################
# Loading data (converted) from MCMax3D 
system_spectrum=np.loadtxt(sim+'/spectrum_PDS70_system.dat')
stellar_spectrum=np.loadtxt(sim+'/spectrum_PDS70_star.dat')

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


############################################################
# Loading profiles modeled
odata=np.loadtxt("alma_radial_profile_observed.dat")
r_obs=odata[:,0:1]
b_obs=odata[:,1:2]
db_obs=odata[:,2:3]

odata_j_1=np.loadtxt("jband_radial_cut_0FWHM_smoothed.dat")
r_obs_j_1=odata_j_1[:,0:1]
b_obs_j_1=odata_j_1[:,1:2]


mprofile_alma=np.loadtxt(sim+"/alma_radial_profile_modeled.dat")
r_alma=np.reshape(mprofile_alma[:,0:1],mprofile_alma.shape[0])
b_alma=np.reshape(mprofile_alma[:,1:2],mprofile_alma.shape[0])

mprofile_jband=np.loadtxt(sim+"/jband_radial_profile_modeled.dat")
r_jband_complete=np.reshape(mprofile_jband[:,0:1],mprofile_jband.shape[0])
b_jband_complete=np.reshape(mprofile_jband[:,1:2],mprofile_jband.shape[0])

r_jband=r_jband_complete
b_jband=b_jband_complete


############################################################
# Working out the Qphi profile
r_jband=[]
b_jband=[]
for i in range(0,len(r_jband_complete)):
    if abs(r_jband_complete[i])>=5.6:
        r_jband.append(r_jband_complete[i])
        b_jband.append(b_jband_complete[i])
r_jband=np.array(r_jband)
b_jband=np.array(b_jband)


############################################################
# Plotting
fig=plt.figure(figsize=(9,3))
gs=gridspec.GridSpec(1,3,wspace=0.0,hspace=0.0)
ax1=plt.subplot(gs[0,0])
ax2=plt.subplot(gs[0,1])
ax3=plt.subplot(gs[0,2])

#fig,(ax1,ax2,ax3)=plt.subplots(1,3,figsize=(8,3))

fsize=12
msize=2.0
lwidth=2.0
mcolor="blue"
ax1.errorbar(r_obs,b_obs,db_obs,fmt='.',label="observation",elinewidth=0.6,markersize=msize,color="black")
ax1.plot(r_alma,b_alma,'-',linewidth=lwidth,label="model",color=mcolor)

ax2.plot(r_obs_j_1,b_obs_j_1,".",label="observation",markersize=msize,color="black")
ax2.plot(r_jband,b_jband,"-",linewidth=lwidth,label="model",color=mcolor)

ax3.plot(x_system,y_system*x_system,linewidth=lwidth,color=mcolor)
ax3.plot(x_star,y_star*x_star,linewidth=lwidth,label="Photospere",color="silver")
ax3.errorbar(x_photometric,y_photometric*x_photometric,y_system_error,fmt='.',capsize=3,color="black")

ax1.set_ylim(0.0,)
ax1.set_xlim(0.0,120.0)
ax2.set_ylim(0.0,4.0)
ax2.set_xlim(-120.0,120.0)
ax3.set_xlim(0.1,)
ax3.set_ylim(1e-17,1e-12)


ax3.set_xscale("log")
ax3.set_yscale("log")

ax1.set_xlabel(r"$r$ (AU)",fontsize=fsize)
ax1.set_ylabel(r"Normalized density flux",fontsize=fsize)
ax2.set_xlabel(r"$r$ (AU)",fontsize=fsize)
ax2.set_ylabel(r"Normalized density flux",fontsize=fsize)
ax3.set_xlabel(r'$\lambda \, (\mu\mathrm{m})$',fontsize=fsize)
ax3.set_ylabel(r'$\lambda F_{\lambda}$ (W/m^2)',fontsize=fsize)

ax1.tick_params(labelsize=fsize)
ax2.tick_params(labelsize=fsize)
ax3.tick_params(labelsize=fsize)

ax1.locator_params(axis='x',nbins=5)
ax2.locator_params(axis='x',nbins=5)

#ax1.legend()
#ax2.legend()
ax3.legend(loc="upper right",fontsize=10)

gs.tight_layout(fig,w_pad=2)
name="rt_model"
#fig.savefig("../%s.png"%(name))
fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/paper_figures/%s.png"%(name))
plt.show()


