
import os
import sys
import math
import mypy
import glob
import time
import pylab
import pickle
import subprocess
from astropy import wcs
from scipy import optimize
from astropy.io import fits
from matplotlib.path import Path
import shapely.geometry as geometry
from scipy.spatial import ConvexHull, Delaunay
from shapely.ops import cascaded_union, polygonize
from scipy.spatial import Voronoi, voronoi_plot_2d
from astropy import wcs, cosmology, constants, units
from matplotlib.backends.backend_pdf import PdfPages


def PolyArea(x, y):
	#  calculate the area of an arbitrary ploygon with given vertices
	return 0.5 * abs(pylab.dot(x, pylab.roll(y, 1)) - pylab.dot(y, pylab.roll(x, 1)))


###  storing basic meta data of EAGLE snapshots
class snapshot:
	def __init__(self, id, expansionFactor, lookbackTime, redshift):
		self.id = id
		self.expansionFactor = expansionFactor
		self.lookbackTime = lookbackTime
		self.redshift = redshift
		self.z = redshift

snapshots_info = [[0, 0.05, 13.62, 20.00],
                 [1, 0.06, 13.53, 15.13],
                 [2, 0.09, 13.32, 9.99],
                 [3, 0.10, 13.25, 8.99],
                 [4, 0.11, 13.16, 8.07],
                 [5, 0.12, 13.04, 7.05],
                 [6, 0.14, 12.86, 5.97],
                 [7, 0.15, 12.75, 5.49],
                 [8, 0.17, 12.63, 5.04],
                 [9, 0.18, 12.46, 4.49],
                 [10, 0.20, 12.25, 3.98],
                 [11, 0.22, 12.01, 3.53],
                 [12, 0.25, 11.66, 3.02],
                 [13, 0.29, 11.16, 2.48],
                 [14, 0.31, 10.86, 2.24],
                 [15, 0.33, 10.53, 2.01],
                 [16, 0.37, 10.05, 1.74],
                 [17, 0.40, 9.49, 1.49],
                 [18, 0.44, 8.86, 1.26],
                 [19, 0.50, 7.93, 1.00],
                 [20, 0.54, 7.37, 0.87],
                 [21, 0.58, 6.71, 0.74],
                 [22, 0.62, 6.01, 0.62],
                 [23, 0.67, 5.19, 0.50],
                 [24, 0.73, 4.16, 0.37],
                 [25, 0.79, 3.23, 0.27],
                 [26, 0.85, 2.29, 0.18],
                 [27, 0.91, 1.34, 0.10],
                 [28, 1.00, 0.00, 0.00]]
snapshots = [snapshot(*info) for info in snapshots_info]

eagle_dir = '../data/eagle_simulations'
sim_catalogs = glob.glob('%s/snap*.dat' % eagle_dir)



###  mass-completeness limits
masslim_catalogs = glob.glob('../data/completeness/masslimits_*_master.dat')
masslim_data = [mypy.readcat(mi) for mi in masslim_catalogs]
masslim95 = pylab.array([mi.masslim95 for mi in masslim_data])
masslim95_zax = masslim_data[0].z
masslim95_average = masslim95.sum(axis=0) / len(masslim95)








###  selecting redshift from EAGLE
si = 20
simdata = mypy.readcat(sim_catalogs[si])


###  reading voronoi map
vmap = fits.getdata('../data/eagle_simulations/medVoronoi.overdens.%02i.fits' % si)
dx_map, dy_map = vmap.shape

overdens_arr = pylab.zeros(len(simdata.GalaxyID))
for i in range(len(simdata.GalaxyID)):
	x, y = int(simdata.x_cMpc[i] * dx_map / 100. - 0.5), int(simdata.y_cMpc[i] * dy_map / 100. - 0.5)
	overdens_arr[i] = vmap[y][x]






###  calculating MFs
dm = 0.25
lmassbins = pylab.arange(9.5-dm/2., 11.5+dm, dm)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.



xlo, xhi = 25. * 100. / dx_map, 100 - 25. * 100. / dx_map
ylo, yhi = 25. * 100. / dy_map, 100 - 25. * 100. / dy_map

inds0 = pylab.find((overdens_arr > 0.0) & (overdens_arr < 0.5) & 
	               (simdata.x_cMpc > xlo) & (simdata.x_cMpc < xhi) & (simdata.y_cMpc > ylo) & (simdata.y_cMpc < yhi))
inds1 = pylab.find((overdens_arr > 0.5) & (overdens_arr < 1.0) &
	               (simdata.x_cMpc > xlo) & (simdata.x_cMpc < xhi) & (simdata.y_cMpc > ylo) & (simdata.y_cMpc < yhi))
inds2 = pylab.find((overdens_arr > 1.0) & (overdens_arr < 1.5) &
	               (simdata.x_cMpc > xlo) & (simdata.x_cMpc < xhi) & (simdata.y_cMpc > ylo) & (simdata.y_cMpc < yhi))
inds3 = pylab.find((overdens_arr > 1.5) & (overdens_arr < 2.0) &
	               (simdata.x_cMpc > xlo) & (simdata.x_cMpc < xhi) & (simdata.y_cMpc > ylo) & (simdata.y_cMpc < yhi))

digi0 = pylab.digitize(pylab.log10(simdata.stellarMass[inds0]), lmassbins)
digi1 = pylab.digitize(pylab.log10(simdata.stellarMass[inds1]), lmassbins)
digi2 = pylab.digitize(pylab.log10(simdata.stellarMass[inds2]), lmassbins)
digi3 = pylab.digitize(pylab.log10(simdata.stellarMass[inds3]), lmassbins)

ngal0 = pylab.bincount(digi0, minlength=len(lmassbins)+1)[1:-1]
ngal1 = pylab.bincount(digi1, minlength=len(lmassbins)+1)[1:-1]
ngal2 = pylab.bincount(digi2, minlength=len(lmassbins)+1)[1:-1]
ngal3 = pylab.bincount(digi3, minlength=len(lmassbins)+1)[1:-1]



























