import mcmax3dpy.read as mread
import mcmax3dpy.plot as mplot
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib import ticker, patches
import astropy.units as u
import sys

zones=mread.read_zones("../output/")

def hstack_matrix(zones,fieldname,vlim=[None,None]):
  k=0
  T=np.zeros((zones[0].np,len(zones)*zones[0].nr))
  R=np.zeros((zones[0].np,len(zones)*zones[0].nr))
  P=np.zeros((zones[0].np,len(zones)*zones[0].nr))
  for zone in zones:
    field=getattr(zone, fieldname)    
    fieldlog=mplot.plog(field)
    vmin=vlim[0]
    vmax=vlim[1]
    
    levels = ticker.MaxNLocator(nbins=39).tick_values(vmax, vmin)    
    ticks = ticker.MaxNLocator(nbins=6, prune="both").tick_values(vmin, vmax)
  
    """
    Midplane cut
    """  
    x=zone.x[:,int(zone.nt/2),:]
    x=np.append(x,[x[0,:]],axis=0)
    y=zone.y[:,int(zone.nt/2),:]
    y=np.append(y,[y[0,:]],axis=0)
    rr=zone.r[:,int(zone.nt/2),:]
    pp=zone.phi[:,int(zone.nt/2),:]
    # take an average, because it is not exaclty at zero. 
    # FIXME: check if I really have used the correc index
    val=(fieldlog[:,int(zone.nt/2),:]+fieldlog[:,int(zone.nt/2+1),:])/2.0
    #val=np.append(val,[val[0,:]],axis=0)
    T[:,k*60:(k+1)*60]=val
    R[:,k*60:(k+1)*60]=rr
    P[:,k*60:(k+1)*60]=pp
    k+=1
    print(T.shape,R.shape,P.shape)

  plt.imshow(T,vmin=0.1,vmax=4)
  plt.xlabel("Radial grid point")
  plt.ylabel("Azimuthal direction")  
  plt.show()
  return T,R,P


def T_radial(T,R,P,phi,**kwargs):

    class Midgridp:
        def __init__(self,i,j):
            self.i=i
            self.j=j

        def getR(self):
            return R[i][j]
            
        def getPhi(self):
            return P[i][j]

        def getT(self):
            return T[i][j]

    ############################################################
    # Find i index closest to the input phi value
    """
    phi=(phi*u.deg).to(u.rad).value
    phi_values=np.reshape(P[:,0:1],P.shape[0])
    islit=(np.abs(phi_values-phi)).argmin()
    """
    islit=round(np.arctan(yp/xp)/(2*np.pi/zone[0].np))-1
    

    R_cut=np.reshape(R[islit:islit+1,:],R.shape[1])
    T_cut=np.reshape(T[islit:islit+1,:],T.shape[1])        

    if kwargs["smooth"]==False:
        return T_cut,R_cut
    elif kwargs["smooth"]==True:
        # Identify index of grid points falling in the same AU along the slit
        # iff r>10 AU
        comd=[] # Array for common distances
        for j in range(0,len(R_cut)-1):
            if R_cut[j]>=20.0:
                if(int(R_cut[j])==int(R_cut[j+1])):
                    comd_val=int(R_cut[j])
                    print(j,j+1,int(R_cut[j]),int(R_cut[j+1]))
                    if comd_val not in comd:
                        comd.append(comd_val)
        dict_comd={}
        for val in comd:
            dict_comd[val]=[]
            for j in range(0,len(R_cut)):
                if(val==int(R_cut[j])):
                    dict_comd[val].append(j)
        #print(dict_comd)
        lims_dict_comd={}
        for key in dict_comd:
            lims_dict_comd[key]=[]
            min_val=min(dict_comd[key])
            max_val=max(dict_comd[key])            
            lims_dict_comd[key].append(min_val)
            lims_dict_comd[key].append(max_val)
        dict_Tave={}
        for key in lims_dict_comd:
            R_value=key
            min_j=int(lims_dict_comd[R_value][0])
            max_j=int(lims_dict_comd[R_value][1])
            Tave=0
            N=0
            for ii in range(min_j,max_j+1):# Radial average
                Tave+=(np.sum(np.reshape(T[:,ii:ii+1],T.shape[0])))/T.shape[0]
                N+=1
            dict_Tave[R_value]=Tave/N
         
        T_smooth=[]
        R_smooth=[]
        for j in range(0,len(R_cut)):
            if(int(R_cut[j]) in dict_Tave):
                R_smooth.append(int(R_cut[j]))
                T_smooth.append(dict_Tave[int(R_cut[j])])
            else:
                R_smooth.append(R_cut[j])
                T_smooth.append(T_cut[j])
                
        #for j in range(0,len(R_smooth)):
         #   print("%.4f %.4f"%(R_smooth[j],T_smooth[j]))

        return T_smooth,R_smooth

T,R,P=hstack_matrix(zones,"temp",vlim=[1,2000])

T_ave,R_ave=T_radial(T,R,P,33.954434912651394,18.35909447011263,smooth=True)

plt.plot(R_ave,T_ave,".")
plt.xscale("log")
plt.show()

f=open("../tprofile.dat","w")
for i,j in zip(R_ave,T_ave):
  f.write("%.15e %.15e\n"%(i,j))
f.close()
