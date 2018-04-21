
#  This script generates catalogs tabulating values measured
#  from Brian's Voronoi maps of galaxies in the photometric
#  catalogs. These new catalogs are line-matched.

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
from scipy import interpolate
from matplotlib.path import Path
import shapely.geometry as geometry
from scipy.spatial import ConvexHull, Delaunay
from shapely.ops import cascaded_union, polygonize
from astropy import wcs, cosmology, constants, units
from matplotlib.backends.backend_pdf import PdfPages



cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

zgrid = pylab.arange(0.4, 1.4, 0.00001)
vgrid = ((1+zgrid)**2 - 1) / ((1+zgrid)**2 + 1) * constants.c.to(units.km / units.s).value
def vrecess(zi):
	###  return the cosmo recessional velocity at input z
	return pylab.interp(zi, zgrid, vgrid)
def zrecess(vi):
	###  return the cosmo redshift at input v [km/s]
	return pylab.interp(vi, vgrid, zgrid)

def PolyArea(x, y):
	#  calculate the area of an arbitrary ploygon with given vertices
    return 0.5 * abs(pylab.dot(x, pylab.roll(y, 1)) - pylab.dot(y, pylab.roll(x, 1)))



def checkPoint(ref_ra, ref_dec, check_ra, check_dec):
	'''
	This function test the ra, dec from one catalogs in the convex hull from the reference catalog
	ConvexHull reference: 
	http://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html#scipy.spatial.ConvexHull
	example: http://stackoverflow.com/questions/31404658/check-if-points-lies-inside-a-convex-hull

	input:
	  ref_ra, ref_dec ---- ra and dec from reference catalog
	  check_ra, check_dec ---- ra dec from to be checked catalog

	output:
	  return a array:
	  if in, return 1
	  if not, return 0
	'''
	points = pylab.zeros((len(ref_ra), 2))
	points[:,0] = ref_ra
	points[:,1] = ref_dec
	hull = ConvexHull(points)
	hull_path = Path(points[hull.vertices])

	check = pylab.zeros(len(check_ra)) - 1
	for i in range(len(check_ra)):
		check[i] = hull_path.contains_point((check_ra[i], check_dec[i]))
	return check

def alpha_shape(points,alpha):
    """
    Compute the alpha shape (concave hull) of a set
    of points.
 
    @param points: Iterable container of points.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    
    if len(points) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull
 
    def add_edge(edges, edge_points, coords, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(coords[ [i, j] ])
 
    #coords = pylab.array([point.coords[0] for point in points])
    coords = points
 
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
 
        # Lengths of sides of triangle
        a = pylab.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = pylab.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = pylab.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
 
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
 
        # Area of triangle by Heron's formula
        area = pylab.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = a*b*c/(4.0*area)
 
        # Here's the radius filter.
        #print circum_r
        if circum_r < 1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
 
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points

def CHullRandomPoint(poly,npnts):
    ra_LB,dec_LB,ra_UB,dec_UB=poly.bounds
    rand_arr = pylab.random.random(2*npnts)
    ra_rand = (ra_UB-ra_LB)*rand_arr[:npnts] + ra_LB
    dec_rand = (dec_UB-dec_LB)*rand_arr[npnts:] + dec_LB
    testinpoly=pylab.zeros(npnts,dtype='bool')
    for iCR in range(0,npnts):
        testinpoly[iCR]=poly.contains(geometry.Point(ra_rand[iCR],dec_rand[iCR]))
    while pylab.sum(testinpoly)<npnts:
        repinds=True-testinpoly
        rand_arr = pylab.random.random(2*(npnts-pylab.sum(testinpoly)))
        ra_rand[repinds] = (ra_UB-ra_LB)*rand_arr[:npnts-pylab.sum(testinpoly)] + ra_LB
        dec_rand[repinds] = (dec_UB-dec_LB)*rand_arr[npnts-pylab.sum(testinpoly):] + dec_LB
        tmptest=testinpoly[repinds]
        for iCR in range(0,npnts-pylab.sum(testinpoly)):
            tmptest[iCR]=poly.contains(geometry.Point(ra_rand[repinds][iCR],dec_rand[repinds][iCR]))
        testinpoly[repinds]=tmptest
    return ra_rand,dec_rand

def CheckPoints(poly,ras,decs):
    g=[poly.contains(geometry.Point(ras[ip],decs[ip])) for ip in range(0,len(ras))]
    return pylab.arange(len(ras))[pylab.array(g)==True]




def gunzip_read_gzip(path2file, readcat=0, readzout=0, dtype=float, delimiter=None, comments='#'):

	if readcat == readzout:
		raise IOError('Need to specify either readcat or readzout!')

	print '  reading: %s' % os.path.basename(path2file)
	subprocess.call('gunzip %s' % path2file, shell=1)
	if readcat:
		outer = mypy.readcat(path2file[:-3], dtype=dtype, delimiter=delimiter, comments=comments)
	elif readzout:
		outer = mypy.readzout(path2file[:-3])
	subprocess.call('gzip %s' % path2file[:-3], shell=1)

	return outer




voronoi_dir = '/Users/atomczak/DATA/ORELSE/Voronoi'

class massFunction_voronoi_slices:
	def __init__(self, name, zclust, sigmaz, area_zspec_arcmin2, zlos, zhis):
		self.name = name
		self.zclust = zclust
		self.sigmaz = sigmaz
		self.area_zspec_arcmin2 = area_zspec_arcmin2

		self.zlos = zlos                  # lower zlimits for voronoi slices
		self.zhis = zhis                  # upper zlimits for voronoi slices
		self.zbars = (zhis + zlos) / 2.
		self.dz = (zhis - zlos) / 2.

		self.dm = 0.25
		self.lmassbins = pylab.arange(8.0-self.dm/2., 11.5+self.dm, self.dm)
		self.lmassbars = (self.lmassbins[1:] + self.lmassbins[:-1]) / 2.

		self.areas_full_arcmin2 = pylab.zeros(len(self.zbars))               # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_field_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_00v05_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_05v10_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_10v15_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_15v20_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice

		self.volumes_full_Mpc3 = pylab.zeros(len(self.zbars))               # volumes corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_field_Mpc3 = pylab.zeros(len(self.zbars))              # volumes corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_00v05_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_05v10_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_10v15_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_15v20_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice

		self.ngals_full = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_full = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_full = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_full = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_field = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_field = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_field = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_field = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.d_overdens = 0.5
		self.overdens_bins = pylab.array([0., 0.5, 1., 1.5, 2.])
		self.overdens_bars = (self.overdens_bins[1:] + self.overdens_bins[:-1]) / 2.
		self.ngals_00v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_00v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_00v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_00v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_05v10 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_05v10 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_05v10 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_05v10 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_10v15 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_10v15 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_10v15 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_10v15 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_15v20 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_15v20 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_15v20 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_15v20 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.areas_all_arcmin2 = [self.areas_00v05_arcmin2, self.areas_05v10_arcmin2, self.areas_10v15_arcmin2, self.areas_15v20_arcmin2]
		self.volumes_all_Mpc3 = [self.volumes_00v05_Mpc3, self.volumes_05v10_Mpc3, self.volumes_10v15_Mpc3, self.volumes_15v20_Mpc3]
		self.ngals_all = [self.ngals_00v05, self.ngals_05v10, self.ngals_10v15, self.ngals_15v20]
		self.elo_ngals_all = [self.elo_ngals_00v05, self.elo_ngals_05v10, self.elo_ngals_10v15, self.elo_ngals_15v20]
		self.ehi_ngals_all = [self.ehi_ngals_00v05, self.ehi_ngals_05v10, self.ehi_ngals_10v15, self.ehi_ngals_15v20]
		self.phis_all = [self.phis_00v05, self.phis_05v10, self.phis_10v15, self.phis_15v20]


images = ['/Users/atomczak/GoogleDrive/ORELSE/images/nep200_v0.0.5/nep200_v0.0.5_detection_ri.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/sc1324_v0.0.2/sc1324_v0.0.2_detection_i.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/rcs0224-0002_v0.0.2/rcs0224-0002_v0.0.2_detection_I+.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/rxj1716+6708_v0.0.7/rxj1716+6708_v0.0.7_detection_RcI+Z+.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/nep5281_v0.0.2/nep5281_v0.0.2_detection_Y.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/sc1604_v0.0.3/sc1604_v0.0.3_detection_Rc.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/cl0910+5422_v0.0.3/cl0910+5422_v0.0.3_detection_RcI+Z+.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/sc0849+4452_v0.0.2/sc0849+4452_v0.0.2_detection_Z+.fits']






class voronoi_catalog:
	def __init__(self, version, ngals, nzbins):
		self.version = version            # photometric catalog
		self.zlos = pylab.zeros(nzbins)
		self.zhis = pylab.zeros(nzbins)
		self.zbars = pylab.zeros(nzbins)
		self.ntot = pylab.zeros(nzbins)   # Median number of galaxies in bin from all iterat
		self.zspec_fraction = pylab.zeros(nzbins)   # Median fraction of galaxies in bin with spec-zs
		self.dens_matrix = pylab.zeros((ngals, nzbins)) - 99
		self.overdens_matrix = pylab.zeros((ngals, nzbins)) - 99
		self.bunit_dens = ''
		self.bunit_overdens = ''

class substructure:
	def __init__(self, name, ra, dec, zmean, vsigma05mpc, err_vsigma05mpc, Ngal05mpc, vsigma1mpc, err_vsigma1mpc, Ngal1mpc, lmvir, err_lmvir):
		'''
		All information comes from ORELSE.dispersions
		name                    name of substructure
		ra                      luminosity weighted central RA
		dec                     luminosity weighted central Dec
		zmean                   mean redshift of substructure
		vsigma05mpc      km/s   velocity dispersion for galaxies at <0.5 Mpc
		err_vsigma05mpc  km/s   error on previous
		Ngal05mpc               number of zspec galaxies at <0.5 Mpc
		vsigma1mpc       km/s   velocity dispersion for galaxies at <1.0 Mpc
		err_vsigma1mpc   km/s   error on previous
		Ngal1mpc                number of zspec galaxies at <1.0 Mpc
		lmvir            lmsol  log of M_vir (see Lemaux+2012)
		err_lmvir        lmsol  error on previous
		'''
		self.name = name
		self.ra = ra
		self.dec = dec
		self.zmean = zmean
		self.vsigma05mpc = vsigma05mpc
		self.vsigma1mpc = vsigma1mpc
		self.err_vsigma05mpc = err_vsigma05mpc
		self.err_vsigma1mpc = err_vsigma1mpc
		self.Ngal05mpc = Ngal05mpc
		self.Ngal1mpc = Ngal1mpc
		self.lmvir = lmvir
		self.err_lmvir = err_lmvir

		self.vmean = vrecess(zmean)
		self.vlo_3sig = self.vmean - 3*vsigma1mpc
		self.vhi_3sig = self.vmean + 3*vsigma1mpc
		self.zlo_3sig = zrecess(self.vlo_3sig)
		self.zhi_3sig = zrecess(self.vhi_3sig)

github_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'
class field:
	def __init__(self, name, version, zclust, sigmaz, alpha=50, chi2red_thresh=10):
		self.name = name          # e.g. "NEP 200"
		self.version = version    # e.g. "nep200_v0.0.4"
		self.zclust = zclust      # cluster redshift
		self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
		self.zspec_lo = 0         # lower redshift bound for specz
		self.zspec_hi = 0         # upper redshift bound for specz
		self.substructures = []
		self.alpha = alpha

		self.cat = gunzip_read_gzip('%s/%s/%s.cat.gz' % (github_dir, version, version), readcat=1)
		self.zout = gunzip_read_gzip('%s/%s/%s.zout.gz' % (github_dir, version, version), readzout=1)
		self.fout = gunzip_read_gzip('%s/%s/%s.fout.gz' % (github_dir, version, version), readcat=1)
		self.restframe = gunzip_read_gzip('%s/%s/%s.restframe.gz' % (github_dir, version, version), readcat=1)
		self.rfcolors = gunzip_read_gzip('%s/%s/%s.restframe_colors.gz' % (github_dir, version, version), readcat=1)


		###  SETTING OBJECTS IDENTIFIED AS SECURE STARS FROM SPECTROSCOPY TO use=0
		self.crossmatch = gunzip_read_gzip('%s/%s/%s.crossmatch.gz' % (github_dir, version, version), readcat=1, dtype=str)
		self.star_inds = pylab.find(self.crossmatch.Q == '-1')
		for i_star in self.star_inds:
			id_phot_arr = self.crossmatch.id_phot[i_star].split(',')
			for id_phot in id_phot_arr:
				if id_phot != '-1':
					self.cat.use[int(id_phot)-1] *= 0


		print '  reading: %s_voronoi.pickle' % name
		self.voronoi = pickle.load(open('../data/%s_voronoi.pickle' % name, 'rb'))
		self.overdens_max = []
		for vi in range(len(self.voronoi.overdens_matrix)):
			self.overdens_max.append(self.voronoi.overdens_matrix[vi].max())
		self.overdens_max = pylab.array(self.overdens_max)


		###  UPDATING USE FLAG WITH REDUCE CHI**2 THRESHOLD
		chi2red = self.zout.chi_p / (self.zout.nfilt - 1.)
		cinds = pylab.find((chi2red > chi2red_thresh) & (self.cat.z_spec < 0))
		self.cat.use[cinds] = 0.


		###  UVJ classification
		self.UV_color = -2.5 * pylab.log10(self.restframe.restflux_U / self.restframe.restflux_V)
		self.VJ_color = -2.5 * pylab.log10(self.restframe.restflux_V / self.restframe.restflux_J)
		self.uvj_class = mypy.uvj_select(self.UV_color, self.VJ_color, self.fout.z)

		###  UVJ classification based on EAZY in 2-filter mode
		self.UV_color_2filter = self.rfcolors.U_V_color
		self.VJ_color_2filter = self.rfcolors.V_J_color
		self.uvj_class_2filter = mypy.uvj_select(self.UV_color_2filter, self.VJ_color_2filter, self.fout.z)

		###  NUVrJ classification
		self.NUVr_color = -2.5 * pylab.log10(self.restframe.restflux_NUV / self.restframe.restflux_r)
		self.rJ_color = -2.5 * pylab.log10(self.restframe.restflux_r / self.restframe.restflux_J)
		self.nuvrj_class = mypy.nuvrj_select(self.NUVr_color, self.rJ_color, self.fout.z)



		xyrd1 = self.cat.x[0], self.cat.y[0], self.cat.ra[0], self.cat.dec[0]
		xyrd2 = self.cat.x[1], self.cat.y[1], self.cat.ra[1], self.cat.dec[1]
		d_arcsec = mypy.radec_sep(xyrd1[2], xyrd1[3], xyrd2[2], xyrd2[3])
		d_pixel = ((xyrd1[0]-xyrd2[0])**2 + (xyrd1[1] - xyrd2[1])**2)**0.5 
		self.px_scale = d_arcsec / d_pixel

		### getting z band magnitude
		try: self.zmagnitude = 25 - 2.5 * pylab.log10(self.cat.fluxauto_z)
		except: pass
		if name != 'SC1324':
			try: self.zmagnitude = 25 - 2.5 * pylab.log10(self.cat.fluxauto_Zplus)
			except: pass

		###  setting spatial flag based on Convex Hull method
		#zspecinds = pylab.find(self.cat.z_spec > 0)
		#self.spatial_flag = checkPoint(self.cat.ra[zspecinds], self.cat.dec[zspecinds], self.cat.ra, self.cat.dec)
		#self.inds_spatial = pylab.find(self.spatial_flag == 1)

		###  setting spatial flag based on Concave Hull method (creates a "tighter" spatial selection than Cnvex Hull)
		zspecinds = pylab.find(self.cat.z_spec > 0)
		points = pylab.array(zip(self.cat.x[zspecinds], self.cat.y[zspecinds]))
		self.concave_hull, self.edge_points = alpha_shape(points, alpha)
		self.inds_spatial = CheckPoints(self.concave_hull.buffer(10), self.cat.x, self.cat.y)
		self.area_pix2 = self.concave_hull.buffer(10).area
		self.area_arcmin2 = self.area_pix2 * (self.px_scale/60.)**2
		print ''


fields = []
fields.append(field('N200',    'nep200_v0.0.5',       0.691,  0.027, alpha=1./600, chi2red_thresh=7))
fields.append(field('SC1324',  'sc1324_v0.0.2',       0.755,  0.033, alpha=1./600, chi2red_thresh=10))
fields.append(field('RCS0224', 'rcs0224-0002_v0.0.2', 0.772,  0.027, alpha=1./500, chi2red_thresh=10))
fields.append(field('RXJ1716', 'rxj1716+6708_v0.0.7', 0.813,  0.021, alpha=1./500, chi2red_thresh=8))
fields.append(field('N5281',   'nep5281_v0.0.2',      0.818,  0.029, alpha=1./1000, chi2red_thresh=10))
#fields.append(field('SG0023',  'sg0023+0423_v0.1.9',  0.845,  0.025, alpha=, chi2red_thresh=14))
fields.append(field('SC1604',  'sc1604_v0.0.3',       0.910,  0.029, alpha=1./500, chi2red_thresh=10))
fields.append(field('SC0910',  'cl0910+5422_v0.0.3',  1.110,  0.035, alpha=1./500, chi2red_thresh=10))
fields.append(field('SC0849',  'sc0849+4452_v0.0.2',  1.261,  0.029, alpha=1./600, chi2red_thresh=10))
print ''



for i in [12203, 13820]: fields[0].cat.use[i-1] *= 0
for i in [4386, 10388, 10907]: fields[1].cat.use[i-1] *= 0
for i in [17351, 26627, 38501, 43517, 51926, 55037, 55089, 61118]: fields[2].cat.use[i-1] *= 0
for i in []: fields[3].cat.use[i-1] *= 0
for i in [11331, 17375]: fields[4].cat.use[i-1] *= 0
for i in [105231, 170586]: fields[5].cat.use[i-1] *= 0
for i in [51114, 57254]: fields[6].cat.use[i-1] *= 0
for i in [9419, 15619, 18725, 22467, 25342, 28122]: fields[7].cat.use[i-1] *= 0



fields[0].substructures.append(substructure('RXJ1757', 269.33196, 66.525991, 0.6931, 541.1, 139.2, 17, 862.3, 107.9, 34, 14.832, 0.250))

fields[1].substructures.append(substructure('Cluster-A*', 201.20097, 30.1920, 0.7556, 1019.2, 142.0, 25, 873.4, 110.8, 43, 14.833, 0.254))
fields[1].substructures.append(substructure('Cluster-B', 201.08830, 30.2156, 0.6979, 897.0, 394.7, 11, 677.1, 143.6, 13, 14.516, 0.424))
fields[1].substructures.append(substructure('Group-C', 201.00770, 30.4164, 0.7382, 143.8, 40.7, 6, 205.9, 90.1, 8, 12.955, 0.875))
fields[1].substructures.append(substructure('Cluster-I', 201.20096, 30.9680, 0.6957, 646.1, 113.0, 16, 891.5, 128.9, 27, 14.875, 0.291))

###  NOTE: The 3rd substructure is a serendipitous
###        group/cluster behind the main LSS
fields[2].substructures.append(substructure('RCS0224-A*+', 36.15714, -0.0949, 0.7780, 713.3, 179.3, 14, 825.4, 193.2, 34, 14.754, 0.468))
fields[2].substructures.append(substructure('RCS0224-B', 36.14123, -0.0394, 0.7781, 815.5, 138.5, 24, 710.7, 58.8, 52, 14.559, 0.165))
#fields[2].substructures.append(substructure('RCS0224-S**', 36.32021, -0.0928, 0.8454, 204.2, 107.7, 7, 437.7, 115.8, 15, 13.911, 0.529))

fields[3].substructures.append(substructure('RXJ1716-A*', 259.20162, 67.1392, 0.8116, 1150.1, 161.8, 25, 1150.2, 113.4, 62, 15.178, 0.197))
fields[3].substructures.append(substructure('RXJ1716-B+', 259.25229, 67.1533, 0.8127, 767.1, 112.0, 20, 693.5, 61.0, 40, 14.519, 0.176))

fields[4].substructures.append(substructure('RXJ1821', 275.38451, 68.465768, 0.8168, 1146.1, 124.8, 27, 1119.6, 99.6, 52, 15.142, 0.178))

fields[5].substructures.append(substructure('Cluster-A', 241.09311, 43.0821, 0.8984, 576.9, 120.8, 23, 722.4, 134.5, 35, 14.551, 0.372))
fields[5].substructures.append(substructure('Cluster-B', 241.10796, 43.2397, 0.8648, 792.9, 80.1, 28, 818.4, 74.2, 49, 14.722, 0.269))
fields[5].substructures.append(substructure('Group-C*', 241.03142, 43.2679, 0.9344, 1079.6, 289.6, 12, 453.5, 39.6, 32, 13.935, 0.039))
fields[5].substructures.append(substructure('Cluster-D+', 241.14094, 43.3539, 0.9227, 675.6, 180.1, 40, 688.2, 88.1, 70, 14.481, 0.256))
fields[5].substructures.append(substructure('Group-F++', 241.20104, 43.3684, 0.9331, 619.9, 135.0, 14, 541.9, 110.0, 20, 14.168, 0.406))
fields[5].substructures.append(substructure('Group-G', 240.92745, 43.4030, 0.9019, 398.5, 85.4, 7, 539.3, 124.0, 18, 14.169, 0.460))
fields[5].substructures.append(substructure('Group-H', 240.89890, 43.3669, 0.8528, 283.2, 71.7, 9, 287.0, 68.3, 10, 13.359, 0.074))
fields[5].substructures.append(substructure('Group-I', 240.79746, 43.3915, 0.9024, 163.0, 65.1, 5, 333.0, 129.4, 7, 13.541, 0.777))

fields[6].substructures.append(substructure('0910-A', 137.51190, 54.3103, 1.1024, 506.8, 235.9, 10, 893.0, 285.8, 18, 14.777, 0.640))
fields[6].substructures.append(substructure('0910-B', 137.68463, 54.3736, 1.1016, 906.3, 192.4, 11, 795.5, 138.1, 22, 14.627, 0.347))

fields[7].substructures.append(substructure('SC0849A+', 132.23463, 44.761780, 1.2622, 609.1,  178.1, 8,  708.0, 186.9, 13, 14.437, 0.344))
fields[7].substructures.append(substructure('SC0849B+', 132.29977, 44.865903, 1.2636, 433.3,  389.7, 6,  286.3, 64.9,  9,  13.257, 0.295))
fields[7].substructures.append(substructure('SC0849C+', 132.24443, 44.866012, 1.2631, 1006.8, 310.8, 17, 766.2, 97.6,  25, 14.540, 0.166))	
fields[7].substructures.append(substructure('SC0849D',  132.15110, 44.899400, 1.2700, 930.9,  147.7, 15, 788.7, 136.1, 24, 14.576, 0.225))
fields[7].substructures.append(substructure('SC0849E+', 132.27496, 44.959253, 1.2600, 366.0,  205.9, 10, 414.1, 107.8, 14, 13.739, 0.339))

###  Note: I'm removing sg0023 because we feel that it isn't belong in this 
###        study due, in large part, to the fact that it's a lower-mass system.
#fields[].substructures.append(substructure('0023-A*', 6.02560, 4.3590, 0.8396, 507.0, 125.5, 10, 412.8, 119.2, 14, 13.836, 0.578))
#fields[].substructures.append(substructure('0023-B1**', 5.97570, 4.3884, 0.8290, 106.2, 51.4, 12, 176.3, 29.6, 23, 12.730, 0.336))
#fields[].substructures.append(substructure('0023-B2+', 5.96970, 4.3820, 0.8453, 231.3, 53.8, 18, 277.8, 41.0, 38, 13.319, 0.295))
#fields[].substructures.append(substructure('0023-C', 5.92470, 4.3807, 0.8466, 543.8, 58.7, 22, 385.3, 54.3, 45, 13.744, 0.282))
#fields[].substructures.append(substructure('0023-M++', 5.96740, 4.3199, 0.8472, 487.3, 84.8, 7, 418.8, 68.9, 14, 13.853, 0.330))






t0 = time.time()

dm = 0.25
lmassbins = pylab.arange(9.5-dm/2., 11.75+dm, dm)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.
corr_factor = pylab.array([ 1.167, 1.091, 1.250, 0.917, 1.071, 1.231, 1.000, 1.000, 1.000, 1.000])
corr_factor_sf = pylab.array([ 0.895, 0.876, 0.787, 0.848, 0.756, 0.936, 0.791, 1.067, 1.000, 1.000])
corr_factor_qu = pylab.array([ 0.333, 1.455, 0.929, 0.831, 0.924, 0.892, 0.889, 0.858, 0.979, 0.979])
#corr_factor_sf = pylab.ones(len(corr_factor))
#corr_factor_qu = pylab.ones(len(corr_factor))




for i_mc in range(100):

	###  reading in MC fout files
	print 'reading MC iter %04i fout catalogs ...\n' % i_mc
	for f in fields:
		f.fout = mypy.readcat('../data/mc/mc_fouts/%s/%s_mc%04i.fout' % (f.version, f.version, i_mc))






	##################################
	###  calculating full MFs for LSS
	##################################

	nfinal = pylab.zeros((len(fields), len(lmassbars)))
	nfinal_sf = pylab.zeros((len(fields), len(lmassbars)))
	nfinal_qu = pylab.zeros((len(fields), len(lmassbars)))

	for i in range(len(fields)):

		f = fields[i]
		kpc2arcsec = cosmo.kpc_proper_per_arcmin(f.zclust).value / 60.


		###  finding voronoi slice centered on zclust
		#i_vslice = pylab.find(abs(f.voronoi.zbars - f.zclust) == abs(f.voronoi.zbars - f.zclust).min())[0]
		#vinds = pylab.find(f.voronoi.overdens_matrix[:, i_vslice] < 0.8)
		#f.cat.use[vinds] = 2.



		###  finding lower/upper zphot limits
		zlo_phot = f.zclust - 1.5 * f.sigmaz * (1 + f.zclust)
		zhi_phot = f.zclust + 1.5 * f.sigmaz * (1 + f.zclust)


		###  finding lower/upper zspec limits
		zlo_spec, zhi_spec = f.zclust, f.zclust
		dr_pkpc_arr = pylab.zeros((len(f.substructures), len(f.cat.id)))
		for ssi in range(len(f.substructures)):
			ss = f.substructures[ssi]
			vmean = vrecess(ss.zmean)
			vlo = vmean - 3*ss.vsigma1mpc
			vhi = vmean + 3*ss.vsigma1mpc
			zlo = zrecess(vlo)
			zhi = zrecess(vhi)
			if zlo < zlo_spec: zlo_spec = zlo
			if zhi > zhi_spec: zhi_spec = zhi

			###  calculating projected distance between SS and each galaxy
			dr_arcsec = mypy.radec_sep(ss.ra, ss.dec, f.cat.ra, f.cat.dec)
			dr_pkpc_arr[ssi] += dr_arcsec * kpc2arcsec


		###  identifying projected radius to nearest substructure for all galaxies
		dr_pkpc_min = pylab.zeros(len(f.cat.id))
		for gi in range(len(f.cat.id)):
			dr_pkpc_min[gi] += min(dr_pkpc_arr[:,gi])






		#print '\nUSING EXPANDED ZSPEC WINDOW !!!\n'
		#zlo_spec = f.zclust - 0.2
		#zhi_spec = f.zclust + 0.2






		###################
		###  ALL GALAXIES
		###################

		###  LSS members selected by zspec
		members_zspec = pylab.find((f.cat.use[f.inds_spatial] == 1) &
	    		                   (dr_pkpc_min[f.inds_spatial] <= 1500.) &
			                       (f.cat.z_spec[f.inds_spatial] > zlo_spec) &
	    		                   (f.cat.z_spec[f.inds_spatial] < zhi_spec))

		###  LSS candidates selected by zphot
		members_zphot = pylab.find((f.cat.use[f.inds_spatial] == 1) &
			                       (f.cat.z_spec[f.inds_spatial] < 0) &
	    		                   (dr_pkpc_min[f.inds_spatial] <= 1500.) &
			                       (f.fout.z[f.inds_spatial] > zlo_phot) &
	   			                   (f.fout.z[f.inds_spatial] < zhi_phot))


		###  binning galaxies by stellar mass
		digi_mass_zspec = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zspec], lmassbins)
		digi_mass_zphot = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zphot], lmassbins)

		ngal_bins_zspec = pylab.bincount(digi_mass_zspec, minlength=len(lmassbins)+1)[1:-1]
		ngal_bins_zphot = pylab.bincount(digi_mass_zphot, minlength=len(lmassbins)+1)[1:-1]
		ngal_bins_zphot_corr = ngal_bins_zphot * corr_factor
		nfinal[i] += ngal_bins_zspec + ngal_bins_zphot_corr






		###################
		###  STAR-FORMING
		###################

		###  LSS members selected by zspec
		members_zspec_sf = pylab.find((f.cat.use[f.inds_spatial] == 1) &
	    		                   (dr_pkpc_min[f.inds_spatial] <= 1500.) &
			                       (f.cat.z_spec[f.inds_spatial] > zlo_spec) &
	    		                   (f.cat.z_spec[f.inds_spatial] < zhi_spec) &
	    		                   (f.uvj_class[f.inds_spatial] == 1))

		###  LSS candidates selected by zphot
		members_zphot_sf = pylab.find((f.cat.use[f.inds_spatial] == 1) &
			                       (f.cat.z_spec[f.inds_spatial] < 0) &
	    		                   (dr_pkpc_min[f.inds_spatial] <= 1500.) &
			                       (f.fout.z[f.inds_spatial] > zlo_phot) &
	   			                   (f.fout.z[f.inds_spatial] < zhi_phot) &
	    		                   (f.uvj_class[f.inds_spatial] == 1))


		###  binning galaxies by stellar mass
		digi_mass_zspec_sf = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zspec_sf], lmassbins)
		digi_mass_zphot_sf = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zphot_sf], lmassbins)

		ngal_bins_zspec_sf = pylab.bincount(digi_mass_zspec_sf, minlength=len(lmassbins)+1)[1:-1]
		ngal_bins_zphot_sf = pylab.bincount(digi_mass_zphot_sf, minlength=len(lmassbins)+1)[1:-1]
		ngal_bins_zphot_corr_sf = ngal_bins_zphot_sf * corr_factor_sf
		nfinal_sf[i] += ngal_bins_zspec_sf + ngal_bins_zphot_corr_sf






		################
		###  QUIESCENT
		################

		###  LSS members selected by zspec
		members_zspec_qu = pylab.find((f.cat.use[f.inds_spatial] == 1) &
	    		                   (dr_pkpc_min[f.inds_spatial] <= 1500.) &
			                       (f.cat.z_spec[f.inds_spatial] > zlo_spec) &
	    		                   (f.cat.z_spec[f.inds_spatial] < zhi_spec) &
	    		                   (f.uvj_class[f.inds_spatial] == 0))

		###  LSS candidates selected by zphot
		members_zphot_qu = pylab.find((f.cat.use[f.inds_spatial] == 1) &
			                       (f.cat.z_spec[f.inds_spatial] < 0) &
	    		                   (dr_pkpc_min[f.inds_spatial] <= 1500.) &
			                       (f.fout.z[f.inds_spatial] > zlo_phot) &
	   			                   (f.fout.z[f.inds_spatial] < zhi_phot) &
	    		                   (f.uvj_class[f.inds_spatial] == 0))


		###  binning galaxies by stellar mass
		digi_mass_zspec_qu = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zspec_qu], lmassbins)
		digi_mass_zphot_qu = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zphot_qu], lmassbins)

		ngal_bins_zspec_qu = pylab.bincount(digi_mass_zspec_qu, minlength=len(lmassbins)+1)[1:-1]
		ngal_bins_zphot_qu = pylab.bincount(digi_mass_zphot_qu, minlength=len(lmassbins)+1)[1:-1]
		ngal_bins_zphot_corr_qu = ngal_bins_zphot_qu * corr_factor_qu
		nfinal_qu[i] += ngal_bins_zspec_qu + ngal_bins_zphot_corr_qu








		###  writing output
		outer = open('../data/mc/mc_mfs/MF_%s_LSS_%03i.dat' % (f.name, i_mc), 'w')
		outer.write('#   Ntot = Nzspec + Nzphot_corr\n')
		outer.write('#   Nzspec = number of LSS members selected by zspec\n')
		outer.write('#   Nzphot = number of candidate LSS members selected in broad zphot bin\n')
		outer.write('#   Nzphot_corr = Nzphot corrected by false +/- fractions to estimate the\n')
		outer.write('#                 number of LSS members missed by spectroscopy\n')
		outer.write('# lmass Ntot Nzspec Nzphot Nzphot_corr')
		outer.write('   Nsf Nzspec_sf Nzphot_sf Nzphot_corr_sf')
		outer.write('   Nqu Nzspec_qu Nzphot_qu Nzphot_corr_qu\n')

		for j in range(len(lmassbars)):
			outer.write(' %5.2f' % lmassbars[j])
			outer.write('  %4i' % (ngal_bins_zspec[j] + ngal_bins_zphot_corr[j] + 0.5))
			outer.write('  %4i' % ngal_bins_zspec[j])
			outer.write('  %4i' % ngal_bins_zphot[j])
			outer.write('  %4i' % (ngal_bins_zphot_corr[j] + 0.5))

			outer.write('     ')
			outer.write('  %4i' % (ngal_bins_zspec_sf[j] + ngal_bins_zphot_corr_sf[j] + 0.5))
			outer.write('  %4i' % ngal_bins_zspec_sf[j])
			outer.write('  %4i' % ngal_bins_zphot_sf[j])
			outer.write('  %4i' % (ngal_bins_zphot_corr_sf[j] + 0.5))

			outer.write('     ')
			outer.write('  %4i' % (ngal_bins_zspec_qu[j] + ngal_bins_zphot_corr_qu[j] + 0.5))
			outer.write('  %4i' % ngal_bins_zspec_qu[j])
			outer.write('  %4i' % ngal_bins_zphot_qu[j])
			outer.write('  %4i' % (ngal_bins_zphot_corr_qu[j] + 0.5))
			outer.write('\n')
		outer.close()
		print 'wrote to: ../data/mc/mc_mfs/MF_%s_LSS_mc%04i.dat' % (f.name, i_mc)

	print '\n'

















	################################################
	###  calculating MFs for each Voronoi zslice
	###  NOTE: these MFs do not have a statistical
	###     correction applied. Just number counts
	################################################

	for i in range(len(fields)):

		f = fields[i]
		voronoi_overdens = fits.open('%s/%s.mastermedVoronoi.overdens.100iterations.fits' % (voronoi_dir, f.name))
		mf_voronoi_slices_field = massFunction_voronoi_slices(f.name, f.zclust, f.sigmaz, f.area_arcmin2, f.voronoi.zlos, f.voronoi.zhis)

		detection_wcs = wcs.WCS(images[i])
		detection_image = fits.open(images[i])
		det_map = detection_image[0].data
		xx_det, yy_det = pylab.meshgrid(range(det_map.shape[1]), range(det_map.shape[0]))

		ra1, dec1 = detection_wcs.wcs_pix2world([0, 0], [0, 1], 1)
		px_scale_detection = mypy.radec_sep(ra1[0], dec1[0], ra1[1], dec1[1])

		for j in range(len(voronoi_overdens)):

			mypy.progress_bar(j, len(voronoi_overdens))

			vmap = voronoi_overdens[j].data
			vheader = voronoi_overdens[j].header
			zlo, zhi = vheader['Z1'], vheader['Z2']
			vwcs = wcs.WCS(vheader)
			xx_vor, yy_vor = pylab.meshgrid(range(vmap.shape[1]), range(vmap.shape[0]))
			ra_vor, dec_vor = vwcs.wcs_pix2world(xx_vor, yy_vor, 1)
			ra1, dec1 = vwcs.wcs_pix2world([0, 0], [0, 1], 1)
			px_scale_vor = mypy.radec_sep(ra1[0], dec1[0], ra1[1], dec1[1])


			'''
			interper = interpolate.interp2d(xx_vor, yy_vor, vmap)

			xnew = pylab.linspace(0, vmap.shape[1]-1, vmap.shape[1]*px_scale_vor/px_scale_detection)
			ynew = pylab.linspace(0, vmap.shape[0]-1, vmap.shape[0]*px_scale_vor/px_scale_detection)
			interp_vmap = interper(xnew, ynew)
			'''


			###  performing concave hull on vornoi map
			zspecinds = pylab.find(f.cat.z_spec > 0)

			###  pixel coords of galxies in voronoi map
			x_vor, y_vor = vwcs.wcs_world2pix(f.cat.ra[zspecinds], f.cat.dec[zspecinds], 1)
			points_vor = pylab.array(zip(x_vor, y_vor))

			xpoints, ypoints = [], []
			for x in range(vmap.shape[1]):
				for y in range(vmap.shape[0]):
					xpoints.append(x)
					ypoints.append(y)
			xpoints = pylab.array(xpoints)
			ypoints = pylab.array(ypoints)

			###  creating "segmentation" map for pixels that are
			###  part of the region of spectroscopic coverage
			concave_hull_vor, edge_points_vor = alpha_shape(points_vor, f.alpha * px_scale_vor / px_scale_detection)
			inds_vor = CheckPoints(concave_hull_vor.buffer(1), xpoints, ypoints)
			inds_vor_spatial = (ypoints[inds_vor]-1, xpoints[inds_vor]-1)
			segmap_vor = vmap * 0.
			segmap_vor[inds_vor_spatial] += 1


			###  calculating sky area/volumes of log(overdensity) thresholds of spectroscopic region
			inds_full = pylab.where(segmap_vor == 1)
			inds_field = pylab.where((segmap_vor == 1) & (vmap <= 0.3))
			mf_voronoi_slices_field.areas_full_arcmin2[j] = len(inds_full[0]) * (px_scale_vor/60.)**2
			mf_voronoi_slices_field.areas_field_arcmin2[j] = len(inds_field[0]) * (px_scale_vor/60.)**2
			mf_voronoi_slices_field.volumes_full_Mpc3[j] = mypy.vcomoving_slice(mf_voronoi_slices_field.areas_field_arcmin2[j], zlo, zhi)
			mf_voronoi_slices_field.volumes_field_Mpc3[j] = mypy.vcomoving_slice(mf_voronoi_slices_field.areas_field_arcmin2[j], zlo, zhi)

			vmap_vbin = vmap * 0.
			vmap_vbin[inds_full] = vmap[inds_full]
			#fits.writeto('../data/voronoi_maps/%s_full_z%.4f.fits' % (f.name, mf_voronoi_slices_field.zbars[j]),
			#             vmap_vbin, header=vheader, clobber=1)
			vmap_vbin = vmap * 0.
			vmap_vbin[inds_field] = vmap[inds_field]
			#fits.writeto('../data/voronoi_maps/%s_field_z%.4f.fits' % (f.name, mf_voronoi_slices_field.zbars[j]),
			#             vmap_vbin, header=vheader, clobber=1)

			###  ... for overdensity bins
			for k in range(len(mf_voronoi_slices_field.overdens_bars)):
				inds_vbin = pylab.where((segmap_vor == 1) & 
					                     (vmap >= mf_voronoi_slices_field.overdens_bins[k]) &
					                     (vmap <= mf_voronoi_slices_field.overdens_bins[k+1]))
				mf_voronoi_slices_field.areas_all_arcmin2[k][j] = len(inds_vbin[0]) * (px_scale_vor/60.)**2
				mf_voronoi_slices_field.volumes_all_Mpc3[k][j] = mypy.vcomoving_slice(mf_voronoi_slices_field.areas_all_arcmin2[k][j], zlo, zhi)

				vmap_vbin = vmap * 0.
				vmap_vbin[inds_vbin] = vmap[inds_vbin]
				#fits.writeto('../data/voronoi_maps/%s_%.1fv%.1f_z%.4f.fits' % (f.name, mf_voronoi_slices_field.overdens_bins[k], mf_voronoi_slices_field.overdens_bins[k+1], mf_voronoi_slices_field.zbars[j]),
				#             vmap_vbin, header=vheader, clobber=1)





			###  number counting galaxies in massbins
			inds_full = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
				                   (f.fout.z[f.inds_spatial] >= zlo) &
				                   (f.fout.z[f.inds_spatial] < zhi))
			inds_field = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
				                    (f.fout.z[f.inds_spatial] >= zlo) &
				                    (f.fout.z[f.inds_spatial] < zhi) &
				                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= 0.3))

			digi_mass_full = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_full], mf_voronoi_slices_field.lmassbins)
			ngal_full = pylab.bincount(digi_mass_full, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
			nlo_full, nhi_full = mypy.massfunc.CI(ngal_full)
			mf_voronoi_slices_field.ngals_full[j] += ngal_full
			mf_voronoi_slices_field.elo_ngals_full[j] += nlo_full
			mf_voronoi_slices_field.ehi_ngals_full[j] += nhi_full
			mf_voronoi_slices_field.phis_full[j] += ngal_full / mf_voronoi_slices_field.volumes_full_Mpc3[j] / mf_voronoi_slices_field.dm

			digi_mass_field = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_field], mf_voronoi_slices_field.lmassbins)
			ngal_field = pylab.bincount(digi_mass_field, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
			nlo_field, nhi_field = mypy.massfunc.CI(ngal_field)
			mf_voronoi_slices_field.ngals_field[j] += ngal_field
			mf_voronoi_slices_field.elo_ngals_field[j] += nlo_field
			mf_voronoi_slices_field.ehi_ngals_field[j] += nhi_field
			mf_voronoi_slices_field.phis_field[j] += ngal_field / mf_voronoi_slices_field.volumes_field_Mpc3[j] / mf_voronoi_slices_field.dm

			###  ... for overdensity bins
			for k in range(len(mf_voronoi_slices_field.overdens_bars)):
				inds_vbin =  pylab.find((f.cat.use[f.inds_spatial] == 1) & 
					                    (f.fout.z[f.inds_spatial] >= zlo) &
					                    (f.fout.z[f.inds_spatial] < zhi) &
					                    (f.voronoi.overdens_matrix[f.inds_spatial, j] >= mf_voronoi_slices_field.overdens_bins[k]) &
					                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= mf_voronoi_slices_field.overdens_bins[k+1]))

				###  if no galaxies exist in this overdensity bin do nothing
				if len(inds_vbin) == 0:
					pass
				else:
					digi_mass_vbin = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_vbin], mf_voronoi_slices_field.lmassbins)
					ngal_vbin = pylab.bincount(digi_mass_vbin, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
					nlo_vbin, nhi_vbin = mypy.massfunc.CI(ngal_vbin)
					mf_voronoi_slices_field.ngals_all[k][j] += ngal_vbin
					mf_voronoi_slices_field.elo_ngals_all[k][j] += nlo_vbin
					mf_voronoi_slices_field.ehi_ngals_all[k][j] += nhi_vbin
					mf_voronoi_slices_field.phis_all[k][j] += ngal_vbin / mf_voronoi_slices_field.volumes_all_Mpc3[k][j] / mf_voronoi_slices_field.dm


		pickle.dump(mf_voronoi_slices_field, open('../data/mc/mc_mfs/MF_%s_voronoi_slices_mc%04i.pickle' % (f.name, i_mc), 'wb'))
		print 'wrote to: ../data/mc/mc_mfs/MF_%s_voronoi_slices.pickle' % f.name
	print '\n'























