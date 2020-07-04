import numpy as np
import matplotlib.pyplot as plt
import sys
"""
fig = plt.figure()
x = np.arange(10) / 10
y = (x + 0.1)**2

plt.errorbar(x, y, xerr=0.5, label='xlolims=True')
y = (x + 0.1)**3


plt.legend()
plt.show()

sys.exit()
"""
data=np.loadtxt("../azprofile_alma_obs.dat")
x=np.reshape(data[:,0:1],data.shape[0])
y=np.reshape(data[:,1:2],data.shape[0])
y_err_beam=np.reshape(data[:,2:3],data.shape[0])
y_err_pixel=np.reshape(data[:,3:4],data.shape[0])

fig1=plt.figure()
ax1=plt.axes()
ax1.errorbar(x,y,yerr=y_err_beam,marker=".",fmt="o",capsize=2,elinewidth=0.5)
ax1.set_title("Error per beam")
ax1.set_xlabel("PA (deg)")
ax1.set_ylabel("Intensity")
plt.show()
fig1.savefig("error_per_beam.png")

fig2=plt.figure()
ax2=plt.axes()
ax2.errorbar(x,y,yerr=y_err_pixel,marker=".",fmt="o",capsize=2,elinewidth=0.5)
ax2.set_title("Error per pixel")
ax2.set_xlabel("PA (deg)")
ax2.set_ylabel("Intensity")
plt.show()
fig2.savefig("error_per_pixel.png")
