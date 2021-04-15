import matplotlib.pyplot as plt, numpy as np, sys
from matplotlib import ticker, patches
from astropy import units as u, constants as c
from mcmax3dpy import read as mread, plot as mplot

zones=mread.read_zones("../output_low_mass/")

def plot_cuts_zones(zones,fieldname):
    fig,ax=plt.subplots()
    for zone in zones:
        field=np.log10(getattr(zone,fieldname))
        x=zone.x[0,0:int(zone.nt*0.5),:]
        y=zone.y[0,0:int(zone.nt*0.5),:]
        z=zone.z[0,0:int(zone.nt*0.5),:]
        rho=(x**2+y**2)**0.5
        val=field[0,0:int(zone.nt*0.5),:]
        levels=np.linspace(val.min(),val.max(),100)
        CS=ax.contourf(rho,z/rho,val,levels=levels)
        CL=ax.contour(rho,z/rho,val,levels=[-24,-21,-18,-15,-12],colors="white")
        ax.clabel(CL,fmt='%.1f',fontsize="smaller",inline=False,
                  rightside_up=False,use_clabeltext=True,colors="black")
        fig.colorbar(CS,format="%.1f",label=fieldname)
        ax.set(xlabel="r (AU)",ylabel="z/r",xscale="log",
               xlim=(rho[-1].min(),rho[-1].max()),ylim=(0.00,0.45))
        plt.show()
    return fig

plot_cuts_zones(zones,"rhod")


