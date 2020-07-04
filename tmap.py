import mcmax3dpy.read as mread
import mcmax3dpy.plot as mplot
import numpy
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
from matplotlib import ticker, patches
import astropy.units as u

zones=mread.read_zones("output/")

fig=mplot.plot_cuts_zones(zones, "temp")
plt.show()
