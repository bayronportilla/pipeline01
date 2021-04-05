import mcmax3dpy.read as mread
import mcmax3dpy.plot as mplot
import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib import ticker, patches
import astropy.units as u
import sys

zones=mread.read_zones("../output/")
fig=mplot.plot_midplane_zones_new(zones, "temp",vlim=[1,1600],plotGrid=False)
plt.xlim(25,41)
plt.ylim(9,27)
#fig=mplot.plot_cuts_zones(zones,"temp",vlim=[1,1600],plotGrid=False)
plt.show()
