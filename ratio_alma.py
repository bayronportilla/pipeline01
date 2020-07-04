import numpy as np
import sys


def ratio_SMsm():
    data=np.loadtxt("azprofile_alma_mod.dat")
    r=np.reshape(data[:,0:1],data.shape[0])
    b=np.reshape(data[:,1:2],data.shape[0])
    
    w=50.0
    min1=79.9
    min2=268.0
    max1=153.0
    max2=353.0

    bmin1=[]
    bmin2=[]
    bmax1=[]
    bmax2=[]

    for (i,j) in zip(r,b):
        if (min1-0.5*w)<=i<=(min1+0.5*w):
            bmin1.append(j)
        elif (min2-0.5*w)<=i<=(min2+0.5*w):
            bmin2.append(j)
        elif (max1-0.5*w)<=i<=(max1+0.5*w):
            bmax1.append(j)
        elif (max2-0.5*w)<=i<=(max2+0.5*w):
            bmax2.append(j)

    bmin1_val=min(bmin1)
    bmin2_val=min(bmin2)
    bmax1_val=max(bmax1)
    bmax2_val=max(bmax2)

    f=open("ratios.dat","w")
    f.write("ratio_1=%.5f\n"%(bmax1_val/bmin1_val))
    f.write("ratio_2=%.5f\n"%(bmax2_val/bmin2_val))
    print("ratio_1=%.5f"%(bmax1_val/bmin1_val))
    print("ratio_2=%.5f"%(bmax2_val/bmin2_val))
    
    return None
