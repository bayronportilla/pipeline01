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
    rho0=np.ones((rho.shape[0],rho.shape[1]))*rho[-1]
    z0=np.flip(z[:,-1])
    z00=np.ones((rho0.shape[0],rho0.shape[1]))
    for i in range(0,z00.shape[0]):
        z00[i,:]=z0[i]

    print(z0)
    #print(rho0)
    print(z00/120)
    #sys.exit()
    levels=np.linspace(f.min(),f.max(),100)
    CS=ax.contourf(rho0,z00/120,f,levels=levels,extend='both')
    CL=ax.contour(rho0,z00/120,f,levels=[-23,-20,-17],colors="white")
    CB=fig.colorbar(CS,format="%d",label=r"log $\rho_\mathrm{dust}\, (\mathrm{g/cm}^3)$")
    CB.set_ticks(np.linspace(f.min(),f.max(),5))
    ax.clabel(CL,fmt='%.1f',fontsize="smaller",inline=False,
              rightside_up=False,use_clabeltext=True,colors="black")
    ax.set(xlabel="r (AU)",ylabel="z/r",xscale="linear")#,
           #xlim=(rho[-1].min(),rho[-1].max()),ylim=(0.0,0.6))    
    plt.show()
    return None
plot_cuts_zones(zones,"rhod")


