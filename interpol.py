import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.interpolate import CubicSpline

def modify_density():

    ############################################################
    # Import data
    data_obs=np.loadtxt("alma_radial_profile_observed.dat")
    data_mod=np.loadtxt("alma_radial_profile_modeled.dat")


    ############################################################
    # Creating independent arrays
    r_obs=np.reshape(data_obs[:,0:1],data_obs.shape[0])
    f_obs=np.reshape(data_obs[:,1:2],data_obs.shape[0])
    r_mod=np.reshape(data_mod[:,0:1],data_mod.shape[0])
    f_mod=np.reshape(data_mod[:,1:2],data_mod.shape[0])


    ############################################################
    # Prepare cubic spline interpolation
    cs_obs=CubicSpline(r_obs,f_obs)
    cs_mod=CubicSpline(r_mod,f_mod)


    ############################################################
    # Determine limits for interpolation
    if min(r_obs)<min(r_mod):
        r_min=min(r_obs)
    else:
        r_min=min(r_mod)

    if max(r_obs)>max(r_mod):
        r_max=max(r_obs)
    else:
        r_max=max(r_mod)


    """
    ############################################################
    # Do a check?
    r=np.linspace(r_min,r_max,100)
    plt.plot(r,cs_obs(r),label='observation interpolated',color='blue')
    plt.plot(r,cs_mod(r),label='model interpolated',color='red')
    plt.plot(r_obs,f_obs,'*',color='blue')
    plt.plot(r_mod,f_mod,'*',color='red')
    plt.legend()
    plt.show()
    """


    ############################################################
    # Import initial density file
    data_ini=np.loadtxt("surface_density_PDS70.dat")
    r_ini=np.reshape(data_ini[:,0:1],data_ini.shape[0])
    d_ini=np.reshape(data_ini[:,1:2],data_ini.shape[0])


    ############################################################
    # Built R coefficients
    R_array=[cs_obs(i)/cs_mod(i) for i in r_ini]


    ############################################################
    # Built new density profile
    d_new=[i*j for i,j in zip(d_ini,R_array)]


    """
    ############################################################
    # Plot final result
    plt.plot(r_ini,d_ini,label='initial density')
    plt.plot(r_ini,d_new,label='corrected density')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    """

    file=open("surface_density_PDS70.dat",'w')
    for i in range(0,len(r_ini)):
        file.write("%.15e %.15e\n"%(r_ini[i],d_new[i]))
    file.close()

    return None




