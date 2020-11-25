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

def zone_matrix(zone,Rin,Rout,**kwargs):
  
  # Reading information from 'input' file
  path_input_file='../input.dat'
  infile=open(path_input_file).readlines()
  for line in infile:
    if line.split('=')[0]=='star02:x':
      xp=float(line.split('=')[1])
    elif line.split('=')[0]=='star02:y':
      yp=float(line.split('=')[1])
    elif line.split('=')[0]=='zone4:Rin':
      Rin=float(line.split('=')[1])
    elif line.split('=')[0]=='zone4:Rout':
      Rout=float(line.split('=')[1])
  
  # Distance to planet from star
  D=(xp**2+yp**2)**0.5

  # Distance to CPD's Rout from star
  d=(D**2+Rout**2)**0.5
  
  # theta,theta prime, theta_1 and theta_2
  theta=np.arctan(yp/xp)
  thetap=np.arcsin(Rout/d)
  theta_1=theta-thetap
  theta_2=theta+thetap

  # Find i index closest to the input phi value
  islit=int(round(theta/(2*np.pi/zone.np)))-1
  islit_1=int(round(theta_1/(2*np.pi/zone.np)))-1
  islit_2=int(round(theta_2/(2*np.pi/zone.np)))-1
  Nislit=islit_2-islit_1+1

  print(Nislit)
  #sys.exit()
  """
  # Print zone information
  print("Radial points: %d"%(zone.nr))
  print("theta points: %d"%(zone.nt))
  print("phi points: %d"%(zone.np))

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

  # Return triplet (they all have the same dimension)
  T=val
  R=rr
  P=pp

  # Find i index closest to the input phi value
  islit=int(round(np.arctan(yp/xp)/(2*np.pi/zone.np)))-1

  # Temperature along islit
  R_islit=np.reshape(R[islit:islit+1,:],R.shape[1])
  T_islit=np.reshape(T[islit:islit+1,:],T.shape[1])
    
  plt.imshow(T)
  plt.show()

  plt.plot(R_islit,T_islit,'.')
  plt.show()
  sys.exit()

  if kwargs['smooth']==True:
    std_array=[]
    for i in range(0,T.shape[1]):
      col=np.reshape(T[:,i:i+1],T.shape[0])
      std_array.append(np.std(col))
    plt.plot(std_array,'.')
    plt.show()
    std_cut=float(input("Enter the limit in the Std: "))
    T_ave=0
    R_ave=[]
    k=0
    for i in range(0,T.shape[1]):
      if std_array[i]>=std_cut: 
        T_ave+=np.sum(np.reshape(T[:,i:i+1],T.shape[0]))
        R_ave.append(R[0][i])
        k+=1
    T_ave=T_ave/k
    R_mid=(R_ave[-1]-R_ave[0])*0.5+R_ave[0]
    sys.exit()
    
    delta=1
    T_smooth=0
    N=0
    iarray=[]
    for i in range(0,len(R[0])):
      if R[0][i]<=Rin+delta:
        
        T_smooth+=(np.sum(np.reshape(T[:,i:i+1],T.shape[0]))/T.shape[0])
        N+=1
       
        iarray.append(i)
    R_islit=[]
    T_islit=[]
    for i in range(0,len(R[0])):
      if i not in iarray:
        R_islit.append(i)
        T_islit.append(T[islit:islit+1,i:i+1])
        
    print(T_islit)
    
    plt.plot(R_islit,T_islit)
    plt.show()

    
    sys.exit()
    print(T[:,iarray[0]:iarray[-1]])
    print(T[:,iarray[0]:iarray[-1]].min(),T[:,iarray[0]:iarray[-1]].max())
    sys.exit()
    
    print(T)
    T_smooth=T_smooth/N
    r_smooth=0.5*(R[0][N-1]-R[0][0])+R[0][0]
    
    
    # Print information about smoothening
    print("Interval covered (AU): ",R[0,0:N])
    print("Number of collapsed columns: %d"%N)
    print("Smoothened temperature (K): %.1f"%10**T_smooth)
    #print(r_smooth)
  """
  #return T,R,P
  return None

#zone_matrix(zones[3],33.954434912651394,18.35909447011263,0.007,1.502,smooth=True)

zone_matrix(zones[0],0.05,20.0,smooth=True)
zone_matrix(zones[1],20,41.0,smooth=True)
zone_matrix(zones[2],41,120.0,smooth=True)
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

T_ave,R_ave=T_radial(T,R,P,33.954434912651394,18.35909447011263,smooth=True)

plt.plot(R_ave,T_ave,".")
plt.xscale("log")
plt.show()

f=open("../tprofile.dat","w")
for i,j in zip(R_ave,T_ave):
  f.write("%.15e %.15e\n"%(i,j))
f.close()
