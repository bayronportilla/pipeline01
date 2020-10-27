import numpy as np
import matplotlib.pyplot as plt
import sys
plt.style.use("fancy")


############################################################
# Loading data
obs=np.loadtxt('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/cut_along_c_obs.dat')
mod_inferior=np.loadtxt('/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run06/cut_along_c_mod.dat')
mod_mean=np.loadtxt('/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run03/cut_along_c_mod.dat')
mod_superior=np.loadtxt('/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run07/cut_along_c_mod.dat')


############################################################
# Here it comes "el sable!"
shift=4


############################################################
# Plot
lwidth=2.0
fig=plt.figure()
ax=plt.axes()
ax.axvspan(34.1,43.1,facecolor="olive",alpha=0.15)
ax.axvline(28.6,color="olive",linestyle='--',linewidth=1.0)
ax.axvline(38.6,color="olive",linestyle='-',linewidth=0.5)
ax.axvline(48.6,color="olive",linestyle='--',linewidth=1.0)
ax.text(39,0.17,"PDS 70 c",rotation='vertical',color='olive')
ax.errorbar(obs[:,0:1],obs[:,1:2],obs[:,2:3],fmt='.',label='Observation',color="black")
ax.plot(mod_inferior[:,0:1]-shift,mod_inferior[:,1:2],label=r'$M_\mathrm{CPD}=0.001 \, M_\mathrm{p}$',linewidth=lwidth,color="pink")
ax.plot(mod_mean[:,0:1]-shift,mod_mean[:,1:2],label=r'$M_\mathrm{CPD}=0.01 \, M_\mathrm{p}$',linewidth=lwidth,color="lightcoral")
ax.plot(mod_superior[:,0:1]-shift,mod_superior[:,1:2],label=r'$M_\mathrm{CPD}=0.1 \, M_\mathrm{p}$',linewidth=lwidth,color="brown")
ax.set_xlabel('Heliocentric distance (AU)')
ax.set_ylabel('Flux density (mJy/beam)')
ax.legend(loc="upper left",frameon=False)
ax.set_xlim(20,55)
ax.set_ylim(0.1,0.4)
plt.tight_layout()
name="explore_cpd_mass"
#fig.savefig("../%s.png"%(name))
fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/paper_figures/%s.png"%(name))
plt.show()




