import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import sys

data=np.loadtxt("radial_cut_0FWHM.dat")
x=np.reshape(data[:,0:1],data.shape[0])
y=np.reshape(data[:,1:2],data.shape[0])

xx=[]
yy=[]
for i,j in zip(x,y):
    if i>=5.6 or i<=-5.6:
        xx.append(i)
        yy.append(j)
xx=np.array(xx)
yy=np.array(yy)

porder=3
Ncoeff=21
yyf=savgol_filter(yy,Ncoeff,porder)
vpeak=2.869
plt.plot(xx,yy/vpeak,"+",label="original signal (non-convolved)")
plt.plot(xx,yyf/vpeak,".",label="filtered signal")
plt.xlabel("r [AU]")
plt.ylabel("Normalized flux density")
plt.ylim(-0.5,3.0)
plt.legend(loc="upper right")
plt.show()

f=open("jband_radial_cut_0FWHM_smoothed.dat","w")
for i in range(0,len(xx)):
    f.write("%.5e %.5e\n"%(xx[i],yyf[i]/vpeak))
f.close()
    
    
