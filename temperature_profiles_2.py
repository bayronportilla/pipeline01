import mcmax3dpy.read as mread
import mcmax3dpy.plot as mplot
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib import ticker, patches
import astropy.units as u
import sys
from astropy.stats import median_absolute_deviation
"""
f1, ax1 = plt.subplots()
ax1.plot(range(0,10))
f1.show()
input("Close the figure and press a key to continue")
f2, ax2 = plt.subplots()
ax2.plot(range(10,20))
f2.show()
input("Close the figure and press a key to continue")
#sys.exit()
#x=np.random.normal(0,2,1000)
x=np.random.uniform(size=3)
x=np.append(x,[2.5])
print()

x = y = z = np.arange(0.0,5.0,1.0)
print(x,y,z)
np.savetxt('test.out',(x,y),newline="\n")   # X is an array
"""

#sys.exit()
def plot_statistics(matrix,array,f,visual=None,**kwargs):
  if kwargs["estimator"]=="mean":
    mean=np.mean(array)
    stdv=np.std(array)
    if visual is True:
      fig,((ax1,ax2))=plt.subplots(1,2,figsize=(15,5))
      ax1.imshow(matrix)
      ax1.set_xlabel("radial")
      ax1.set_ylabel("azimuthal")
      ax2.plot(array,".")
      ax2.axhline(mean)
      ax2.axhspan(mean-f*stdv,mean+f*stdv,alpha=0.2)
      ax2.text(0,mean,"Mean=%.1f"%mean)
      ax2.text(0,mean-f*stdv,"$-%d\sigma$"%f)
      ax2.text(0,mean+f*stdv,"$+%d\sigma$"%f)
      ax2.set_xlabel("radial")
      ax2.set_ylabel("standard deviation")
      return fig
    else:
      return mean,stdv
  elif kwargs["estimator"]=="median":
    median=np.median(array)
    mad=median_absolute_deviation(array)
    if visual is True:
      fig,((ax1,ax2))=plt.subplots(1,2,figsize=(15,5))
      ax1.imshow(matrix)
      ax1.set_xlabel("radial")
      ax1.set_ylabel("azimuthal")
      ax2.plot(array,".")
      ax2.axhline(median)
      ax2.axhspan(median-f*mad,median+f*mad,alpha=0.2)
      ax2.text(0,median,"Median=%.1f"%median)
      ax2.text(0,median-f*mad,"$-%d\mathrm{MAD}$"%f)
      ax2.text(0,median+f*mad,"$+%d\mathrm{MAD}$"%f)
      ax2.set_xlabel("radial")
      ax2.set_ylabel("standard deviation")
      return fig
    else:
      return median,mad

  
#plt.show()
    

#fig=plot_statistics(x,1,estimator="median")
#plt.show()
#plot_statistics(x,estimator="mean")
#sys.exit()

zones=mread.read_zones("../output/")

def zone_matrix(zones,ave=False):
  zoneID=0
  for zone in zones:
    zoneID+=1
    # Get temperature field
    field=getattr(zone,"temp")
    fieldlog=mplot.plog(field)

    # Midplane coordinates
    x=zone.x[:,int(zone.nt/2),:]
    y=zone.y[:,int(zone.nt/2),:]
    x=np.append(x,[x[0,:]],axis=0)
    y=np.append(y,[y[0,:]],axis=0)
    rr=zone.r[:,int(zone.nt/2),:]
    pp=zone.phi[:,int(zone.nt/2),:]

    # Midplane temperature (averaged)
    val=(fieldlog[:,int(zone.nt/2),:]+fieldlog[:,int(zone.nt/2+1),:])/2.0

    # Return triplet 
    T,R,P=val,np.reshape(rr[0],rr.shape[0]),np.reshape(pp[:,0:1],pp.shape[0])
    #plt.imshow(T)
    #plt.show()

    if ave is True:
      # Compute array of standard deviation
      std_array=[]
      for i in range(0,T.shape[1]):
        col=np.reshape(T[:,i:i+1],T.shape[0])
        std_array.append(np.std(col))  
        #std_array.append(np.median(col))
      while True:
        f=int(input("Enter the f-value: "))
        fig=plot_statistics(T,std_array,f,estimator="median",visual=True)
        #fig=plot_statistics(T,std_array,f,estimator="mean",visual=True)
        fig.show()
        usef=bool(int(input("Use this f? (1 means yes, 0 means No): ")))
        if usef==True:
          break
      print("Ok, I will use f=",f)
      median=plot_statistics(T,std_array,f,estimator="median",visual=False)[0]
      mad=plot_statistics(T,std_array,f,estimator="median",visual=False)[1]
      T_ave,r_ave=[],[]
      file=open("temp_ave_zone%d.dat"%zoneID,"w")
      for i in range(0,T.shape[1]):
        if abs(std_array[i]-median)<f*mad and std_array[i]!=0:
          T_ave.append(np.sum(np.reshape(T[:,i:i+1],T.shape[0]))/T.shape[0])
          r_ave.append(R[i])
      T_ave,r_ave=np.array(T_ave),np.array(r_ave)
      #np.savetxt("temp_ave_zone%d.dat"%zoneID,(r_ave,T_ave))
      for i,j in zip(r_ave,T_ave):
        file.write("%.15e %.15e\n"%(i,j))
      file.close()
      fig2,ax2=plt.subplots()
      ax2.plot(r_ave,10**T_ave,".")    
      ax2.set_xlabel("r (AU)")
      ax2.set_ylabel("T (K)")
      ax2.set_xscale("log")
      ax2.set_yscale("log")
      plt.show()
      
      print("line")
      #sys.exit()
        
    if ave is False:
      file=open("temp_nonave_zone%d.dat"%zoneID,"w")
      T_cut=np.reshape(T[int(T.shape[0]*0.5):int(T.shape[0]*0.5+1)],T.shape[0])
      for i,j in zip(R,T_cut):
        file.write("%.15e %.15e\n"%(i,j))
      file.close()
      
      

  sys.exit()
    
  if zoneID!=4:

    # Find temperature submatrix of azimuthal directions
    T_sub=T[islit_1:islit_2+1,:]
    #T_sub=T

    # Initialize relevant variables
    T_noisy=0
    j_noisy=[]
    T_clean=[]
    R_clean=[]
    k=0

    
    # Find imax
    rmax=Rins[zoneID-1]+1.0
    for i in range(0,T.shape[1]):
      if rmax>=R[0][i]:
        imax=i
    

    for i in range(0,T.shape[1]):
      selector=abs(std_array[i]-median)/mad
      if selector>=3.0 and i<=imax:
        T_noisy+=np.sum(np.reshape(T[:,i:i+1],T.shape[0]))/T.shape[0]
        j_noisy.append(i)
      elif selector<3.0:
        T_clean.append(np.sum(np.reshape(T_sub[:,i:i+1],T_sub.shape[0]))/T_sub.shape[0])
        R_clean.append(R[0][i])

    if len(j_noisy)!=0:
      T_noisy=T_noisy/len(j_noisy)
      R_noisy=(R[0][j_noisy[-1]]-R[0][j_noisy[0]])*0.5+R[0][j_noisy[0]]
      T_clean.append(T_noisy)
      R_clean.append(R_noisy)


    T_tot=np.array(T_clean)
    R_tot=np.array(R_clean)

    return R_tot,T_tot

  else:
    T_cpd=[]
    R_cpd=[]
    R_cpd_sym=[]
    for i in range(0,T.shape[1]):
      T_cpd.append(np.sum(np.reshape(T[:,i:i+1],T.shape[0]))/T.shape[0])
      R_cpd.append(R[0][i])

    plt.plot(R_cpd,T_cpd,'.')
    plt.show()

    f=open("temp_cpd.dat","w")
    for i in range(0,len(R_cpd)):
      f.write("%.15e %.15e\n"%(R_cpd[i],T_cpd[i]))
    f.close()

    return R_cpd,T_cpd

  return None


zone_matrix(zones,ave=True)
sys.exit()
R_def=[]
T_def=[]
for i in range(1,4):
  R=zone_matrix(i)[0]
  T=zone_matrix(i)[1]
  for j in range(0,len(R)):
    R_def.append(R[j])
    T_def.append(T[j])
  

f=open("temp.dat","w")
for i in range(0,len(T_def)):
  f.write("%.15f %.15f\n"%(R_def[i],T_def[i]))
f.close()

  

sys.exit()


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


def T_radial(T,R,P,xp,yp,**kwargs):

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
    islit=int(round(np.arctan(yp/xp)/(2*np.pi/zones[0].np)))-1


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

T_noisy,R_ave=T_radial(T,R,P,33.954434912651394,18.35909447011263,smooth=True)

plt.plot(R_ave,T_noisy,".")
plt.xscale("log")
plt.show()

f=open("../tprofile.dat","w")
for i,j in zip(R_ave,T_noisy):
  f.write("%.15e %.15e\n"%(i,j))
f.close()
