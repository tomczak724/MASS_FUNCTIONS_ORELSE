
###  This script is designed to track the time-evolution
###  of the Voronoi overdensity of galaxies in the 
###  EAGLE simulation.

import time
import mypy
import pylab
from astropy.io import fits
from astropy import cosmology

cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

###  storing basic meta data of EAGLE snapshots
class snapshot:
	def __init__(self, id, expansionFactor, lookbackTime, redshift):
		self.id = id
		self.expansionFactor = expansionFactor
		self.lookbackTime = lookbackTime
		self.redshift = redshift
		self.z = redshift

snapshots_info = [[15, 0.33, 10.53, 2.01],
                  [16, 0.37, 10.05, 1.74],
                  [17, 0.40, 9.49, 1.49],
                  [18, 0.44, 8.86, 1.26],
                  [19, 0.50, 7.93, 1.00],
                  [20, 0.54, 7.37, 0.87],
                  [21, 0.58, 6.71, 0.74],
                  [22, 0.62, 6.01, 0.62],
                  [23, 0.67, 5.19, 0.50],
                  [24, 0.73, 4.16, 0.37]]
snapshots = [snapshot(*info) for info in snapshots_info]
redshifts = pylab.array([snp.redshift for snp in snapshots])
age = cosmo.age(redshifts).value



###  reading catalogs
eagle_dir = '../data/eagle_simulations'
cats = [mypy.readcat('%s/snap%02i.dat' % (eagle_dir, snp.id)) for snp in snapshots]

###  reading overdensity maps
vmaps = [fits.getdata('%s/medVoronoi.overdens.%02i.fits' % (eagle_dir, snp.id)) for snp in snapshots]







###  indices of galaxies to track
inds0 = pylab.find((cats[0].stellarMass))
pylab.shuffle(inds0)







x_all = pylab.zeros((len(inds0), len(snapshots)))
y_all = pylab.zeros((len(inds0), len(snapshots)))

x_min = pylab.zeros(len(inds0))
y_min = pylab.zeros(len(inds0))
x_max = pylab.zeros(len(inds0))
y_max = pylab.zeros(len(inds0))

lmass_all = pylab.zeros((len(inds0), len(snapshots)))
overdens_all = pylab.zeros((len(inds0), len(snapshots)))




t0 = time.time()

for i, ind0 in enumerate(inds0):

	mypy.progress_bar(i, len(inds0))

	###  tracking galaxy through the first N-1 snapshots
	ind_tracker = ind0 * 1
	for j in range(len(cats)-1):
		cat_b4 = cats[j]
		cat_afta = cats[j+1]

		x_all[i][j] = cat_b4.x_cMpc[ind_tracker]
		y_all[i][j] = cat_b4.y_cMpc[ind_tracker]
		lmass_all[i][j] = pylab.log10(cat_b4.stellarMass[ind_tracker])

		vmap = vmaps[j]
		lx, ly = vmap.shape
		xtran = int(cat_b4.x_cMpc[ind_tracker] * lx / 100. - 0.5)   # pixel coordinate of galaxy in voronoi map
		ytran = int(cat_b4.y_cMpc[ind_tracker] * ly / 100. - 0.5)   # pixel coordinate of galaxy in voronoi map
		overdens_all[i][j] = vmap[ytran][xtran]

		try:
			ind_tracker = pylab.find(cat_afta.GalaxyID == cat_b4.DescendantID[ind_tracker])[0]
		except:
			break


	###  adding final snapshot
	x_all[i][j+1] = cat_b4.x_cMpc[ind_tracker]
	y_all[i][j+1] = cat_b4.y_cMpc[ind_tracker]
	lmass_all[i][j+1] = pylab.log10(cat_b4.stellarMass[ind_tracker])

	vmap = vmaps[j+1]
	lx, ly = vmap.shape
	xtran = int(cat_b4.x_cMpc[ind_tracker] * lx / 100. - 0.5)   # pixel coordinate of galaxy in voronoi map
	ytran = int(cat_b4.y_cMpc[ind_tracker] * ly / 100. - 0.5)   # pixel coordinate of galaxy in voronoi map
	overdens_all[i][j+1] = vmap[ytran][xtran]


	###  grabbing min/max spatial coordinates of each galaxy
	x_min[i] = x_all[i].min()
	y_min[i] = y_all[i].min()
	x_max[i] = x_all[i].max()
	y_max[i] = y_all[i].max()


tf = time.time()








