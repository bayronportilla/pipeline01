import numpy as np
import matplotlib.pyplot as plt
import sys

def find_ratios():


    # Loading data
    data=np.loadtxt("jband_radial_profile_modeled.dat")
    x=np.reshape(data[:,0:1],data.shape[0])
    y=np.reshape(data[:,1:2],data.shape[0])
    

    # Exctracting positive radius
    x_plus=[]
    y_plus=[]
    for i,j in zip(x,y):
        if i>=0:
            x_plus.append(i)
            y_plus.append(j)


    # Finding nearest points
    x_plus=np.array(x_plus)
    y_plus=np.array(y_plus)
    A_obs=(18.5,2.61)
    B_obs=(38.0,0.51)
    C_obs=(54.6,1.00)    
    A_mod=(x_plus[np.abs(x_plus-A_obs[0]).argmin()],y_plus[np.abs(x_plus-A_obs[0]).argmin()])
    B_mod=(x_plus[np.abs(x_plus-B_obs[0]).argmin()],y_plus[np.abs(x_plus-B_obs[0]).argmin()])
    C_mod=(x_plus[np.abs(x_plus-C_obs[0]).argmin()],y_plus[np.abs(x_plus-C_obs[0]).argmin()])
    

    # Printing info
    print("\nPrinting error in relative local minima's heights")

    D1_obs=2.1
    D2_obs=0.49
    D1_mod=A_mod[1]-B_mod[1]
    D2_mod=C_mod[1]-B_mod[1]
    
    print("D1_error=%.2f percent"%(abs(D1_obs-D1_mod)/D1_obs * 100.0))
    print("D2_error=%.2f percent"%(abs(D2_obs-D2_mod)/D2_obs * 100.0))
    print()

    f=open("ratios_jband_radial_flux.dat","w")
    f.write("D1_error=%.2f percent \n"%(abs(D1_obs-D1_mod)/D1_obs * 100.0))
    f.write("D2_error=%.2f percent \n"%(abs(D2_obs-D2_mod)/D2_obs * 100.0))
    f.close()

    return None
