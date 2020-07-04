import numpy as np 
import matplotlib.pyplot as plt
import sys

z1_in=0.04
z2_in=20.0
z3_in=50.0
z1_out=z2_in
z2_out=z3_in
z3_out=120.0


Npoints=1000
r_array=np.linspace(z1_in,z3_out,Npoints)

params_a={'sh_z1':18.8,'sh_z2':15.9,'sh_z3':10.4,
         'shpow_z1':1.15,'shpow_z2':1.15,'shpow_z3':1.15,
         'rsh_z1':120,'rsh_z2':120,'rsh_z3':120}

params_b={'sh_z1':19.0,'sh_z2':25.62,'sh_z3':14.0,
         'shpow_z1':1.15,'shpow_z2':1.7,'shpow_z3':1.5,
          'rsh_z1':120,'rsh_z2':120,'rsh_z3':120}




def sh(r,params):
    if z1_in<=r<z1_out:
        value=params['sh_z1']*(r/params['rsh_z1'])**params['shpow_z1']
    elif z2_in<=r<z2_out:
        value=params['sh_z2']*(r/params['rsh_z2'])**params['shpow_z2']
    else:
        value=params['sh_z3']*(r/params['rsh_z3'])**params['shpow_z3']
    return value

H_array_plus_a=np.array([sh(i,params_a) for i in r_array])
H_array_minus_a=-H_array_plus_a

H_array_plus_b=np.array([sh(i,params_b) for i in r_array])
H_array_minus_b=-H_array_plus_b

fig=plt.figure()


plt.plot(r_array,H_array_plus_a,color='grey',label='1B')
plt.plot(r_array,H_array_minus_a,color='grey')
plt.plot(r_array,H_array_plus_b,color='orange',label='paper')
plt.plot(r_array,H_array_minus_b,color='orange')
plt.axvline(z1_out,linestyle='--',color='blue',linewidth=0.3)
plt.axvline(z2_out,linestyle='--',color='blue',linewidth=0.3)
plt.axvline(z3_out,linestyle='--',color='blue',linewidth=0.3)
plt.legend()
plt.show()


    
