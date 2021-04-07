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
    field=getattr(zone,"temp") # [phi,theta,r]
    fieldlog=mplot.plog(field)

    # Midplane coordinates
    x=zone.x[:,int(zone.nt/2),:]
    y=zone.y[:,int(zone.nt/2),:]
    x=np.append(x,[x[0,:]],axis=0)
    y=np.append(y,[y[0,:]],axis=0)
    rr=zone.r[:,int(zone.nt/2),:] 
    pp=zone.phi[:,int(zone.nt/2),:]

    # Midplane temperature (averaged)
    T=(fieldlog[:,int(zone.nt/2),:]+fieldlog[:,int(zone.nt/2+1),:])*0.5
    R,P=np.reshape(rr[0],rr.shape[1]),np.reshape(pp[:,0:1],pp.shape[0])
    plt.imshow(T)
    plt.show()
    print(zoneID)

    if ave is True:
      # Compute array of standard deviation
      std_array=[]
      for i in range(0,T.shape[1]):
        col=np.reshape(T[:,i:i+1],T.shape[0])
        std_array.append(np.std(col))  
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
      #T_cut=np.reshape(T[int(T.shape[0]*0.5):int(T.shape[0]*0.5+1)],T.shape[1])
      T_cut=np.reshape(T[int(15):int(15+1)],T.shape[1])
      for i,j in zip(R,T_cut):
        file.write("%.15e %.15e\n"%(i,j))
      file.close()
      
    #return None


zone_matrix(zones,ave=False)



