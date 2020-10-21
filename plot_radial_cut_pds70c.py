import numpy as np
import matplotlib.pyplot as plt
import sys
plt.style.use("fancy")


############################################################
# Loading data
obs=np.loadtxt('/Users/users/bportilla/Documents/first_project/scripts/PDS70/observations/cut_along_c_obs.dat')
mod05=np.loadtxt('/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run02/cut_along_c_mod.dat')
mod5=np.loadtxt('/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run03/cut_along_c_mod.dat')
mod12=np.loadtxt('/data/users/bportilla/runs/final_runs/peregrine/1BS_CPD_mass_corrected/run04/cut_along_c_mod.dat')


############################################################
# Here it comes "el sable!"
shift=4


############################################################
# Plot
lwidth=2.0
fig=plt.figure()
ax=plt.axes()
ax.axvspan(34.1,43.1,facecolor="olive",alpha=0.15)
ax.axvline(38.6,color="olive",linestyle='-',linewidth=0.5)
ax.text(35.8,0.85,"PDS 70 c",rotation='vertical',color='olive')
ax.errorbar(obs[:,0:1],obs[:,1:2],obs[:,2:3],fmt='.',label='Observation',color="black")
ax.plot(mod05[:,0:1]-shift,mod05[:,1:2],label=r'$0.5 \, M_{\mathrm{Jup}}$',linewidth=lwidth,color="darkgrey")
ax.plot(mod5[:,0:1]-shift,mod5[:,1:2],label=r'$5 \, M_{\mathrm{Jup}}$',linewidth=lwidth,color="lightcoral")
ax.plot(mod12[:,0:1]-shift,mod12[:,1:2],label=r'$12 \, M_{\mathrm{Jup}}$',linewidth=lwidth,color="peru")
ax.set_xlabel('Heliocentric distance (AU)')
ax.set_ylabel('Flux density (mJy/beam)')
ax.legend(loc="upper left",frameon=False)
ax.set_xlim(0,120)
ax.set_ylim(0,)
plt.tight_layout()
name="cpd_effect_three_cases"
#fig.savefig("../%s.png"%(name))
fig.savefig("/Users/users/bportilla/Documents/first_project/scripts/PDS70/paper_figures/%s.png"%(name))
plt.show()




