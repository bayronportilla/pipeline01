import matplotlib.pyplot as plt, numpy as np, sys
from matplotlib import ticker, patches
from astropy import units as u, constants as c
from mcmax3dpy import read as mread, plot as mplot
plt.style.use("fancy")

zones=mread.read_zones("../output/")

def plot_cuts_zones(zones,fieldname):
    x=np.empty((int(zones[0].nt*0.5),int(zones[0].nr)),dtype=object)
    y,z,f=x,x,x
    fig,ax=plt.subplots()
    for zone in zones:
        field=np.log10(getattr(zone,fieldname))
        x=np.append(x,zone.x[0,0:int(zone.nt*0.5),:],axis=1)
        y=np.append(y,zone.y[0,0:int(zone.nt*0.5),:],axis=1)
        z=np.append(z,zone.z[0,0:int(zone.nt*0.5),:],axis=1)
        f=np.append(f,field[0,0:int(zone.nt*0.5),:],axis=1)
    x,y,z=x[:,int(zones[0].nr):],y[:,int(zones[0].nr):],z[:,int(zones[0].nr):]
    f=f[:,int(zones[0].nr):]
    rho=(x**2+y**2)**0.5
    levels=np.linspace(f.min(),f.max(),100)
    CS=ax.contourf(rho,z/rho,f,levels=levels)
    CL=ax.contour(rho,z/rho,f,levels=[-24,-21,-18,-15,-12],colors="white")
    ax.clabel(CL,fmt='%.1f',fontsize="smaller",inline=False,
              rightside_up=False,use_clabeltext=True,colors="black")
    fig.colorbar(CS,format="%.1f",label=fieldname)
    ax.set(xlabel="r (AU)",ylabel="z/r",xscale="log",
           xlim=(rho[-1].min(),rho[-1].max()),ylim=(0.0,0.6))    
    plt.show()
    return None#fig
plot_cuts_zones(zones,"rhod")


