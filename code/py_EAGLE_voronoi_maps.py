
import glob
import time
import mypy
import pylab
import matplotlib
from astropy import wcs
from scipy import signal, ndimage
from astropy.io import fits
from matplotlib import collections
import matplotlib.patches as patches
from scipy.spatial import Voronoi, voronoi_plot_2d

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





###  mass-completeness limits
masslim_catalogs = glob.glob('../data/completeness/masslimits_*_master.dat')
masslim_data = [mypy.readcat(mi) for mi in masslim_catalogs]
masslim95 = pylab.array([mi.masslim95 for mi in masslim_data])
masslim95_zax = masslim_data[0].z
masslim95_average = masslim95.sum(axis=0) / len(masslim95)








eagle_dir = '../data/eagle_simulations'
sim_catalogs = glob.glob('%s/snap*.dat' % eagle_dir)

t0 = time.time()
for si in range(15, 24+1):

	simdata = mypy.readcat(sim_catalogs[si])

	###  selecting galaxies above the average ORELSE mass limit
	masslim = pylab.interp(snapshots[si].redshift, masslim95_zax, masslim95_average)
	inds = pylab.find(simdata.stellarMass > 10**masslim)



	###  generating Voronoi cells
	points = pylab.array([[simdata.x_cMpc[i], simdata.y_cMpc[i]] for i in inds])
	vor = Voronoi(points)



	'''
	###  plotting voronoi grid
	for reg in vor.regions:

		if len(reg) > 0:
			reg.append(reg[0])
			inds = pylab.find(pylab.array(reg) > -1)
			for i in range(len(inds) - 1):

				xy1, xy2 = vor.vertices[reg[inds[i]]], vor.vertices[reg[inds[i+1]]]
				pl.plot([xy1[0], xy2[0]], [xy1[1], xy2[1]], lw=1, color='gray')
	'''



	###  calculating areas of Voronoi cells
	areas = []
	for ri in vor.point_region:
		region = vor.regions[ri]
		x_polygon = pylab.array([vor.vertices[i][0] for i in region])
		y_polygon = pylab.array([vor.vertices[i][1] for i in region])
		areas.append(PolyArea(x_polygon, y_polygon))
	areas = pylab.array(areas)




	###  initial Voronoi map with 75x75 ckpc pixel scale
	#dx, dy = int(100 * (1000 / 75.)), int(100 * (1000 / 75.))
	#xgrid, ygrid = pylab.meshgrid(pylab.linspace(0, 100, dx), pylab.linspace(0, 100, dy))
	#vmap0 = pylab.zeros((dy, dx))

	#print '\n\ngenerating map for %s' % sim_catalogs[si]
	#for xi in range(dx):
	#	mypy.progress_bar(xi, dx)
	#	for yi in range(dy):
	#		dr = ((vor.points[:,0] - xgrid[yi][xi])**2 + (vor.points[:,1] - ygrid[yi][xi])**2)**0.5
	#		i_rmin = pylab.find(dr == dr.min())[0]
	#		vmap0[yi][xi] = 1. / areas[i_rmin]

	#overdens_map = pylab.log10(vmap0 / pylab.median(vmap0))
	#fits.writeto('../data/eagle_simulations/medVoronoi.overdens.%02i.fits' % si, overdens_map, clobber=1)

tf = time.time()

print 'time elapsed: %.1f mins' % ((tf-t0)/60.)











