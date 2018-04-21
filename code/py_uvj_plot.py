
#  This script generates catalogs tabulating values measured
#  from Brian's Voronoi maps of galaxies in the photometric
#  catalogs. These new catalogs are line-matched.

import os
import sys
import math
import mypy
import glob
import time
import ezgal
import pylab
import numpy
import pickle
import matplotlib
import subprocess
from astropy import wcs
from scipy import optimize
from astropy.io import fits
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




def uvj_select(uv, vj, z, uvj_slope=0.88, uvj_intercept=0.59):
    """
    #  INPUT:
    #    uv  --->  Restframe (U-V) color (in AB system)  [array-like]
    #    vj  --->  Restframe (V-J) color (in AB system)  [array-like]
    #    z  ---->  Galaxy redshifts  [array-like]
    #
    #  OUTPUT: quiescent=0, star-forming=1
    #
    #  EXAMPLE:
    #            [ SF, SF, Qui, SF,  ....  Qui, Qui]
    #     array( [ 1,  1,  0,   1,   ....  0,   0 ] )
    #
    #  based on Whitaker+2011
    """
    if type(uv)==list or type(uv)==float or type(uv)==int or type(uv)==long or type(uv)==numpy.float64: uv = numpy.array(uv)
    if type(vj)==list or type(vj)==float or type(vj)==int or type(vj)==long or type(vj)==numpy.float64: vj = numpy.array(vj)

    floor = pylab.zeros(len(z))
    wall = pylab.zeros(len(z))
    slope = pylab.zeros(len(z))
    intercept = pylab.zeros(len(z))

    floor[ pylab.where( (0.0<=z) & (z<1.5) ) ] = 1.3
    floor[ pylab.where( (1.5<=z) & (z<2.0) ) ] = 1.3
    floor[ pylab.where( (2.0<=z) ) ] = 1.2

    wall[ pylab.where( (0.0<=z) & (z<1.5) ) ] = 1.6
    wall[ pylab.where( (1.5<=z) & (z<2.0) ) ] = 1.5
    wall[ pylab.where( (2.0<=z) ) ] = 1.4
    wall += 10

    slope[ pylab.where( z<0.5 ) ] = uvj_slope
    slope[ pylab.where( 0.5<=z ) ] = uvj_slope

    intercept[ pylab.where( z<0.5 ) ] = uvj_intercept + 0.1
    intercept[ pylab.where( 0.5<=z ) ] = uvj_intercept

    outer = pylab.zeros(len(z))
    outer[ pylab.where( (uv<slope*vj+intercept) | (uv<floor) | (vj>wall) )[0] ] = 1
    return outer.astype(int)




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





def schechter_mf(xaxis, alpha, mstar, phistar):
    """
    #  DESCRIPTION:
    #    Returns the values for a Schechter mass function
    #    from a given mass-axis and Schechter parameters.
    #
    #  INPUTS:
    #    xaxis = input mass value(s)
    #    alpha = Schechter parameters
    #    mstar = Schechter parameters
    #    phistar = Schechter parameters    
    """
    return pylab.log(10) * phistar * 10**((xaxis-mstar)*(1+alpha)) * pylab.exp(-10**(xaxis-mstar))

def dschechter(lmassax, lmstar, a1, a2, phistar1, phistar2):
	factor1 = pylab.log(10) * pylab.exp(-10**(lmassax - lmstar)) * 10**(lmassax - lmstar)
	factor2 = phistar1 * 10**(a1*(lmassax - lmstar)) + phistar2 * 10**(a2*(lmassax - lmstar))
	return factor1 * factor2

def smf_TOT_parameterized_leja2015(lmassax, zi):
	'''
	Returns the number densities at the input M* values
	based on the parameterized SMFs from Leja+2015. See
	their Section 2.
	'''
	a1, a2 = -0.39, -1.53
	lphistar1 = -2.64 + 0.07*zi - 0.28*zi**2
	lphistar2 = -3.11 - 0.18*zi - 0.03*zi**2
	lmstar = 10.72 - 0.13*zi + 0.11*zi**2
	return dschechter(lmassax, lmstar, a1, a2, 10**lphistar1, 10**lphistar2)





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
	def __init__(self, name, version, zclust, sigmaz, alpha=50, chi2red_thresh=10, uvj_slope=0.88, uvj_intercept=0.59):
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
		self.fout = gunzip_read_gzip('%s/%s/%s_ZFparams.fout.gz' % (github_dir, version, version), readcat=1)
		self.restframe = gunzip_read_gzip('%s/%s/%s.restframe.gz' % (github_dir, version, version), readcat=1)
		#self.rfcolors = gunzip_read_gzip('%s/%s/%s.restframe_colors.gz' % (github_dir, version, version), readcat=1)


		###  SETTING OBJECTS IDENTIFIED AS SECURE STARS FROM SPECTROSCOPY TO use=0
		self.crossmatch = gunzip_read_gzip('%s/%s/%s.crossmatch.gz' % (github_dir, version, version), readcat=1, dtype=str)
		self.star_inds = pylab.find(self.crossmatch.Q == '-1')
		for i_star in self.star_inds:
			id_phot_arr = self.crossmatch.id_phot[i_star].split(',')
			for id_phot in id_phot_arr:
				if id_phot != '-1':
					self.cat.use[int(id_phot)-1] *= 0


		print '  reading: %s_voronoi.pickle' % name
		self.voronoi = pickle.load(open('/Users/atomczak/Dropbox/ORELSE/DATA/VORONOI/%s_voronoi.pickle' % name, 'rb'))
		self.overdens_max = []
		for vi in range(len(self.voronoi.overdens_matrix)):
			self.overdens_max.append(self.voronoi.overdens_matrix[vi].max())
		self.overdens_max = pylab.array(self.overdens_max)

		###  mass completeness limits
		self.masslimits = mypy.readcat('../data/completeness/masslimits_%s_master.dat' % name)

		###  UPDATING USE FLAG WITH REDUCE CHI**2 THRESHOLD
		chi2red = self.zout.chi_p / (self.zout.nfilt - 1.)
		cinds = pylab.find((chi2red > chi2red_thresh) & (self.cat.z_spec < 0))
		self.cat.use[cinds] = 0.


		###  UVJ classification
		self.UV_color = -2.5 * pylab.log10(self.restframe.restflux_U / self.restframe.restflux_V)
		self.VJ_color = -2.5 * pylab.log10(self.restframe.restflux_V / self.restframe.restflux_J)

		self.UV_color_err = pylab.sqrt((2.5 / pylab.log(10))**2 * \
			                           ((self.restframe.errflux_U / self.restframe.restflux_U)**2 + \
			                           	(self.restframe.errflux_V / self.restframe.restflux_V)**2))
		self.VJ_color_err = pylab.sqrt((2.5 / pylab.log(10))**2 * \
			                           ((self.restframe.errflux_V / self.restframe.restflux_V)**2 + \
			                           	(self.restframe.errflux_J / self.restframe.restflux_J)**2))

		self.uvj_class = uvj_select(self.UV_color, self.VJ_color, self.fout.z, uvj_slope=uvj_slope, uvj_intercept=uvj_intercept)

		###  UVJ classification based on EAZY in 2-filter mode
		#self.UV_color_2filter = self.rfcolors.U_V_color
		#self.VJ_color_2filter = self.rfcolors.V_J_color
		#self.uvj_class_2filter = mypy.uvj_select(self.UV_color_2filter, self.VJ_color_2filter, self.fout.z)

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
fields.append(field('N200',    'nep200_v0.0.7',       0.691,  0.027, alpha=1./600, chi2red_thresh=7, uvj_intercept=0.49))
fields.append(field('SC1324',  'sc1324_v0.0.4',       0.755,  0.033, alpha=1./600, chi2red_thresh=10))
fields.append(field('RCS0224', 'rcs0224-0002_v0.0.4', 0.772,  0.027, alpha=1./500, chi2red_thresh=10))
fields.append(field('RXJ1716', 'rxj1716+6708_v0.0.8', 0.813,  0.021, alpha=1./500, chi2red_thresh=8))
fields.append(field('N5281',   'nep5281_v0.0.4',      0.818,  0.029, alpha=1./1000, chi2red_thresh=10))
#fields.append(field('SG0023',  'sg0023+0423_v0.2.0',  0.845,  0.025, alpha=, chi2red_thresh=14))
fields.append(field('SC1604',  'sc1604_v0.0.6',       0.910,  0.029, alpha=1./500, chi2red_thresh=10))
#fields.append(field('CL1429',  'cl1429+4241_v0.0.3',  0.920,  0.051, alpha=1./500, chi2red_thresh=10))
fields.append(field('SC0910',  'cl0910+5422_v0.0.5',  1.110,  0.035, alpha=1./500, chi2red_thresh=10))
fields.append(field('SC0849',  'sc0849+4452_v0.0.4',  1.261,  0.029, alpha=1./600, chi2red_thresh=10, uvj_intercept=0.54))
print ''



fields[0].version_old = 'nep200_v0.0.5'
fields[1].version_old = 'sc1324_v0.0.2'
fields[2].version_old = 'rcs0224-0002_v0.0.2'
fields[3].version_old = 'rxj1716+6708_v0.0.5'
fields[3].version_old = 'nep5281_v0.0.2'
fields[4].version_old = 'sc1604_v0.0.4'
fields[5].version_old = 'cl0910+5422_v0.0.3'
fields[6].version_old = 'sc0849+4452_v0.0.2'





'''
for i in [12203, 13820]: fields[0].cat.use[i-1] *= 0
for i in [4386, 10388, 10907]: fields[1].cat.use[i-1] *= 0
for i in [17351, 26627, 38501, 43517, 51926, 55037, 55089, 61118]: fields[2].cat.use[i-1] *= 0
for i in []: fields[3].cat.use[i-1] *= 0
for i in [11331, 17375]: fields[4].cat.use[i-1] *= 0
for i in [105231, 170586]: fields[5].cat.use[i-1] *= 0
for i in [51114, 57254]: fields[6].cat.use[i-1] *= 0
for i in [9419, 15619, 18725, 22467, 25342, 28122]: fields[7].cat.use[i-1] *= 0
'''


ii = 0
fields[ii].substructures.append(substructure('RXJ1757', 269.33196, 66.525991, 0.6931, 541.1, 139.2, 17, 862.3, 107.9, 34, 14.832, 0.250))

ii += 1
fields[ii].substructures.append(substructure('Cluster-A*', 201.20097, 30.1920, 0.7556, 1019.2, 142.0, 25, 873.4, 110.8, 43, 14.833, 0.254))
fields[ii].substructures.append(substructure('Cluster-B', 201.08830, 30.2156, 0.6979, 897.0, 394.7, 11, 677.1, 143.6, 13, 14.516, 0.424))
fields[ii].substructures.append(substructure('Group-C', 201.00770, 30.4164, 0.7382, 143.8, 40.7, 6, 205.9, 90.1, 8, 12.955, 0.875))
fields[ii].substructures.append(substructure('Cluster-I', 201.20096, 30.9680, 0.6957, 646.1, 113.0, 16, 891.5, 128.9, 27, 14.875, 0.291))

###  NOTE: The 3rd substructure is a serendipitous
###        group/cluster behind the main LSS
ii += 1
fields[ii].substructures.append(substructure('RCS0224-A*+', 36.15714, -0.0949, 0.7780, 713.3, 179.3, 14, 825.4, 193.2, 34, 14.754, 0.468))
fields[ii].substructures.append(substructure('RCS0224-B', 36.14123, -0.0394, 0.7781, 815.5, 138.5, 24, 710.7, 58.8, 52, 14.559, 0.165))
#fields[ii].substructures.append(substructure('RCS0224-S**', 36.32021, -0.0928, 0.8454, 204.2, 107.7, 7, 437.7, 115.8, 15, 13.911, 0.529))

ii += 1
fields[ii].substructures.append(substructure('RXJ1716-A*', 259.20162, 67.1392, 0.8116, 1150.1, 161.8, 25, 1150.2, 113.4, 62, 15.178, 0.197))
fields[ii].substructures.append(substructure('RXJ1716-B+', 259.25229, 67.1533, 0.8127, 767.1, 112.0, 20, 693.5, 61.0, 40, 14.519, 0.176))

ii += 1
fields[ii].substructures.append(substructure('RXJ1821', 275.38451, 68.465768, 0.8168, 1146.1, 124.8, 27, 1119.6, 99.6, 52, 15.142, 0.178))

ii += 1
fields[ii].substructures.append(substructure('Cluster-A', 241.09311, 43.0821, 0.8984, 576.9, 120.8, 23, 722.4, 134.5, 35, 14.551, 0.372))
fields[ii].substructures.append(substructure('Cluster-B', 241.10796, 43.2397, 0.8648, 792.9, 80.1, 28, 818.4, 74.2, 49, 14.722, 0.269))
fields[ii].substructures.append(substructure('Group-C*', 241.03142, 43.2679, 0.9344, 1079.6, 289.6, 12, 453.5, 39.6, 32, 13.935, 0.039))
fields[ii].substructures.append(substructure('Cluster-D+', 241.14094, 43.3539, 0.9227, 675.6, 180.1, 40, 688.2, 88.1, 70, 14.481, 0.256))
fields[ii].substructures.append(substructure('Group-F++', 241.20104, 43.3684, 0.9331, 619.9, 135.0, 14, 541.9, 110.0, 20, 14.168, 0.406))
fields[ii].substructures.append(substructure('Group-G', 240.92745, 43.4030, 0.9019, 398.5, 85.4, 7, 539.3, 124.0, 18, 14.169, 0.460))
fields[ii].substructures.append(substructure('Group-H', 240.89890, 43.3669, 0.8528, 283.2, 71.7, 9, 287.0, 68.3, 10, 13.359, 0.074))
fields[ii].substructures.append(substructure('Group-I', 240.79746, 43.3915, 0.9024, 163.0, 65.1, 5, 333.0, 129.4, 7, 13.541, 0.777))

###  CL 1429
#ii += 1


ii += 1
fields[ii].substructures.append(substructure('0910-A', 137.51190, 54.3103, 1.1024, 506.8, 235.9, 10, 893.0, 285.8, 18, 14.777, 0.640))
fields[ii].substructures.append(substructure('0910-B', 137.68463, 54.3736, 1.1016, 906.3, 192.4, 11, 795.5, 138.1, 22, 14.627, 0.347))

ii += 1
fields[ii].substructures.append(substructure('SC0849A+', 132.23463, 44.761780, 1.2622, 609.1,  178.1, 8,  708.0, 186.9, 13, 14.437, 0.344))
fields[ii].substructures.append(substructure('SC0849B+', 132.29977, 44.865903, 1.2636, 433.3,  389.7, 6,  286.3, 64.9,  9,  13.257, 0.295))
fields[ii].substructures.append(substructure('SC0849C+', 132.24443, 44.866012, 1.2631, 1006.8, 310.8, 17, 766.2, 97.6,  25, 14.540, 0.166))	
fields[ii].substructures.append(substructure('SC0849D',  132.15110, 44.899400, 1.2700, 930.9,  147.7, 15, 788.7, 136.1, 24, 14.576, 0.225))
fields[ii].substructures.append(substructure('SC0849E+', 132.27496, 44.959253, 1.2600, 366.0,  205.9, 10, 414.1, 107.8, 14, 13.739, 0.339))

###  Note: I'm removing sg0023 because we feel that it isn't belong in this 
###        study due, in large part, to the fact that it's a lower-mass system.
#fields[].substructures.append(substructure('0023-A*', 6.02560, 4.3590, 0.8396, 507.0, 125.5, 10, 412.8, 119.2, 14, 13.836, 0.578))
#fields[].substructures.append(substructure('0023-B1**', 5.97570, 4.3884, 0.8290, 106.2, 51.4, 12, 176.3, 29.6, 23, 12.730, 0.336))
#fields[].substructures.append(substructure('0023-B2+', 5.96970, 4.3820, 0.8453, 231.3, 53.8, 18, 277.8, 41.0, 38, 13.319, 0.295))
#fields[].substructures.append(substructure('0023-C', 5.92470, 4.3807, 0.8466, 543.8, 58.7, 22, 385.3, 54.3, 45, 13.744, 0.282))
#fields[].substructures.append(substructure('0023-M++', 5.96740, 4.3199, 0.8472, 487.3, 84.8, 7, 418.8, 68.9, 14, 13.853, 0.330))


'''
print '\nreading expanded parameter space FAST catalogs...'
fields[0].fout = mypy.readcat('%s/nep200_v0.0.5/nep200_v0.0.5_ZFparams.fout.gz' % github_dir)
fields[1].fout = mypy.readcat('%s/sc1324_v0.0.2/sc1324_v0.0.2_ZFparams.fout.gz' % github_dir)
fields[2].fout = mypy.readcat('%s/rcs0224-0002_v0.0.2/rcs0224-0002_v0.0.2_ZFparams.fout.gz' % github_dir)
fields[3].fout = mypy.readcat('%s/rxj1716+6708_v0.0.7/rxj1716+6708_v0.0.7_ZFparams.fout.gz' % github_dir)
fields[4].fout = mypy.readcat('%s/nep5281_v0.0.2/nep5281_v0.0.2_ZFparams.fout.gz' % github_dir)
fields[5].fout = mypy.readcat('%s/sc1604_v0.0.3/sc1604_v0.0.3_ZFparams.fout.gz' % github_dir)
#fields[6].fout = mypy.readcat('%s/cl1429+4241_v0.0.1/cl1429+4241_v0.0.1_ZFparams.fout.gz' % github_dir)
fields[6].fout = mypy.readcat('%s/cl0910+5422_v0.0.3/cl0910+5422_v0.0.3_ZFparams.fout.gz' % github_dir)
fields[7].fout = mypy.readcat('%s/sc0849+4452_v0.0.2/sc0849+4452_v0.0.2_ZFparams.fout.gz' % github_dir)
print 'done!\n'
'''



###  estimating change in UVJ colors from Av=1

def calzetti(wavelength, R_v=4.05):

	k_lam = pylab.zeros(len(wavelength))
	x     = wavelength / (1.e4)  # convert wavelength to microns
	x1    = pylab.where(x <  0.63)[0]
	x2    = pylab.where(x >= 0.63)[0]

	k_lam[x1] = 2.659 * (-2.156 + 1.509/x[x1] - 0.198/(x[x1]**2) + 0.011/(x[x1]**3)) + R_v
	k_lam[x2] = 2.659 * (-1.857 + 1.040/x[x2]) + R_v

	return k_lam / R_v



Av = 1.
wavelengths = pylab.arange(2000., 15000., 1)
Alam2Av = calzetti(wavelengths)
corr = 10**(-0.4 * Av * Alam2Av)

gitdir = '/Users/atomczak/GitHub'
res = mypy.sps.read_res('%s/ORELSE/Catalogs/tomczak_catalogs/filters/filters.orelse.v07.res' % gitdir)
u = res[24-1]
v = res[26-1]
j = res[27-1]

dUV_Av1 = -2.5 * pylab.log10(mypy.synphot(pylab.array(zip(wavelengths, corr)), u.data) / mypy.synphot(pylab.array(zip(wavelengths, corr)), v.data))
dVJ_Av1 = -2.5 * pylab.log10(mypy.synphot(pylab.array(zip(wavelengths, corr)), v.data) / mypy.synphot(pylab.array(zip(wavelengths, corr)), j.data))




















t0 = time.time()

fig = pylab.figure(figsize=(13.2, 6.2))
sp1 = fig.add_subplot(121, aspect=1)
sp2 = fig.add_subplot(122, aspect=1)

fig.subplots_adjust(wspace=0.3, left=0.08, right=0.97)

for sp in [sp1, sp2]:
	sp.minorticks_on()

	#sp.set_xlabel('$\mathrm{(} \, V \, - \, J \, \mathrm{)}_{\mathrm{rest}}$', style='italic', size=22)
	#sp.set_ylabel('$\mathrm{(} \, U \, - \, V \, \mathrm{)}_{\mathrm{rest}}$', style='italic', size=22)

	#sp.set_xlabel('$\mathrm{( \, V} \, - \, \mathrm{J \, )}_{\mathrm{rest}}$', size=22)
	#sp.set_ylabel('$\mathrm{( \, U} \, - \, \mathrm{V \, )}_{\mathrm{rest}}$', size=22)

	sp.set_xlabel('(V - J)$_{\mathrm{rest}}$', size=22)
	sp.set_ylabel('(U - V)$_{\mathrm{rest}}$', size=22)

mugshot_names = []
vj_zphot, uv_zphot = [], []
vj_zspec, uv_zspec = [], []

vj_zphot_err, uv_zphot_err = [], []
vj_zspec_err, uv_zspec_err = [], []

for i_f in range(len(fields)):

	f = fields[i_f]

	masslim_galaxies95 = pylab.interp(f.fout.z, f.masslimits.z, f.masslimits.masslim_empirical)
	masslim_galaxies_ssp = pylab.interp(f.fout.z, f.masslimits.z, f.masslimits.masslim_empirical)
	masslim_galaxies = (masslim_galaxies95 + masslim_galaxies_ssp) / 2.

	inds_zphot = pylab.find((f.cat.use == 1) &
		                    (f.cat.z_spec < 0) &
		                    (f.fout.z > 0.55) & 
		                    (f.fout.z < 1.3) &
		                    (f.fout.lmass > masslim_galaxies))

	inds_zspec = pylab.find((f.cat.use == 1) &
		                    (f.cat.z_spec > 0.55) & 
		                    (f.cat.z_spec < 1.3) &
		                    (f.fout.lmass > 0*masslim_galaxies95))


	vj_zphot += f.VJ_color[inds_zphot].tolist()
	uv_zphot += f.UV_color[inds_zphot].tolist()
	vj_zphot_err += f.VJ_color_err[inds_zphot].tolist()
	uv_zphot_err += f.UV_color_err[inds_zphot].tolist()

	vj_zspec += f.VJ_color[inds_zspec].tolist()
	uv_zspec += f.UV_color[inds_zspec].tolist()
	vj_zspec_err += f.VJ_color_err[inds_zspec].tolist()
	uv_zspec_err += f.UV_color_err[inds_zspec].tolist()

	#for ind in inds_zphot:
	#	mugshot_names.append(' ./%s/mugshot%06i_%s.pdf' % (f.name, f.cat.id[ind], f.version_old))



mugshot_names = pylab.array(mugshot_names)
vj_zphot, uv_zphot = pylab.array(vj_zphot), pylab.array(uv_zphot)
vj_zspec, uv_zspec = pylab.array(vj_zspec), pylab.array(uv_zspec)
vj_zphot_err, uv_zphot_err = pylab.array(vj_zphot_err), pylab.array(uv_zphot_err)
vj_zspec_err, uv_zspec_err = pylab.array(vj_zspec_err), pylab.array(uv_zspec_err)





#################################
###  plotting photometric sample
#################################

nbins_zphot = 75
xlims, ylims = [-0.1, 2.1], [0.3, 2.5]
hist2d_zphot, xedges_zphot, yedges_zphot = pylab.histogram2d(vj_zphot, uv_zphot,
	                                       bins=(nbins_zphot, nbins_zphot),
	                                       range=(xlims, ylims))
extent_zphot = [xedges_zphot[0], xedges_zphot[-1], yedges_zphot[0], yedges_zphot[-1]]



nlo, nhi = hist2d_zphot.min(), hist2d_zphot.max()
dn = abs(nlo - nhi)

cmap_colors = [pylab.cm.gray_r(i) for i in pylab.linspace(0, 1, 15)]
cmap = pylab.cm.gray_r
cmap = matplotlib.colors.ListedColormap(cmap_colors, name='custom grays', N=None)

colormap = sp1.imshow(hist2d_zphot.T, extent=extent_zphot, interpolation='nearest', cmap=cmap)
colormap.set_clim(nlo+dn*0.03, nhi-dn*0.03)
colormap.set_clim(nlo+dn*0.02, nhi-dn*0.03)
sp1.axis(xlims + ylims)

#ntext = sp1.text(0.06, 0.86, 'Phot. sample\nN = %i' % len(uv_zphot), transform=sp1.transAxes)
ntext = sp1.text(0.03, 0.96, 'Photometric sample\n0.55 < $z_{\mathrm{phot}}$ < 1.3\nN = %i' % len(uv_zphot),
	             horizontalalignment='left', verticalalignment='top', transform=sp1.transAxes)
ntext.set_text('Photometric sample\n0.55 < $z_{\mathrm{phot}}$ < 1.3\nN = %i' % 87887)






#################################
###  plotting photometric sample
#################################

nbins_zspec = 50
xlims, ylims = [-0.1, 2.1], [0.3, 2.5]
hist2d_zspec, xedges_zspec, yedges_zspec = pylab.histogram2d(vj_zspec, uv_zspec,
	                                       bins=(nbins_zspec, nbins_zspec),
	                                       range=(xlims, ylims))
extent_zspec = [xedges_zspec[0], xedges_zspec[-1], yedges_zspec[0], yedges_zspec[-1]]


nlo, nhi = hist2d_zspec.min(), hist2d_zspec.max()
dn = abs(nlo - nhi)

#cmap_colors = [pylab.cm.gray_r(i) for i in pylab.linspace(0, 1, 15)]
#cmap = matplotlib.colors.ListedColormap(cmap_colors, name='custom grays', N=None)
#cmap = pylab.cm.gray_r


xgrid_zspec, ygrid_zspec = pylab.meshgrid((xedges_zspec[1:] + xedges_zspec[:-1])/2.,
	                                      (yedges_zspec[1:] + yedges_zspec[:-1])/2.)
inds = pylab.where((hist2d_zspec > 0) & 
	               (hist2d_zspec <= 2) &
	               (xgrid_zspec.T < 1.0) &
	               (ygrid_zspec.T > 0.7))
hist2d_zspec[inds] = pylab.randint(1, 4, size=(nbins_zspec, nbins_zspec))[inds]



#colormap = sp2.imshow(hist2d_zspec.T[::-1,:], extent=extent_zspec, interpolation='nearest', cmap=cmap)
colormap = sp2.imshow(hist2d_zspec.T, extent=extent_zspec, interpolation='nearest', cmap=cmap)
colormap.set_clim(nlo+dn*0.03, nhi-dn*0.03)
colormap.set_clim(0, nhi)
sp2.axis(xlims + ylims)

#ntext = sp2.text(0.06, 0.86, 'Spec. sample\nN = %i' % len(uv_zspec), transform=sp2.transAxes)
ntext = sp2.text(0.03, 0.96, 'Spectroscopic sample\n0.55 < $z_{\mathrm{spec}}$ < 1.3\nN = %i' % len(uv_zspec),
	             horizontalalignment='left', verticalalignment='top', transform=sp2.transAxes)







brown = '#86592d'
x0, y0 = 1.25, 0.45
for sp in [sp1, sp2]:

	mypy.uvj_select_region(0.8, subplot=sp, kolor='w', lw=5)
	mypy.uvj_select_region(0.8, subplot=sp, kolor='k', lw=2.5)

	a = sp.annotate('', (x0+dVJ_Av1, y0+dUV_Av1),
                    (x0, y0),
                    ha="right", va="center",
                    size=20,
                    arrowprops=dict(arrowstyle='->',
                    	            lw=2,
                                    shrinkA=5,
                                    shrinkB=5,
                                    fc='k', ec='k'))



	t = sp.text(x0 + dVJ_Av1/2, y0 + dUV_Av1/2., '$\Delta \mathrm{A_V}$=1',
	        horizontalalignment='left', verticalalignment='top')












































##############################################################
##############################################################
###  Plotting median error maps for restframe colors.
###
###  Briefly, these will be 2D histgrams in (V-J) vs. (U-V)
###  phase space, but instead of the colorbar scaling with
###  the NUMBER of objects in each cell it will correspond
###  to the MEDIAN (V-J) error and MEDIAN (U-V) error.
##############################################################
##############################################################


N_map_zphot = hist2d_zphot * 0.
uv_errMap_zphot = hist2d_zphot * 0.
vj_errMap_zphot = hist2d_zphot * 0.


for i_y in range(len(yedges_zphot)-1):

	mypy.progress_bar(i_y, len(yedges_zphot))

	for i_x in range(len(xedges_zphot)-1):

		ylo, yhi = yedges_zphot[i_y], yedges_zphot[i_y+1]
		xlo, xhi = xedges_zphot[i_x], xedges_zphot[i_x+1]

		inds_cell = pylab.find((uv_zphot > ylo) & (uv_zphot < yhi) & \
			                   (vj_zphot > xlo) & (vj_zphot < xhi))

		N_map_zphot[i_y][i_x] = len(inds_cell)
		uv_errMap_zphot[i_y][i_x] = pylab.median(uv_zphot_err[inds_cell])
		vj_errMap_zphot[i_y][i_x] = pylab.median(vj_zphot_err[inds_cell])


		###  setting cell in error maps
		###  to 0 if N_in_cell < 5
		if len(inds_cell) < 5:
			uv_errMap_zphot[i_y][i_x] = 0.
			vj_errMap_zphot[i_y][i_x] = 0.


print ''
















fig = pylab.figure(figsize=(9.3, 15.2))
sp1 = fig.add_subplot(211, aspect=1)
sp2 = fig.add_subplot(212, aspect=1)

fig.subplots_adjust(wspace=0.3, left=0.1, right=0.97, bottom=0.07)

sp1.minorticks_on()
sp1.set_xlabel('(V - J)$_{\mathrm{rest}}$', size=22)
sp1.set_ylabel('(U - V)$_{\mathrm{rest}}$', size=22)
sp1.axis(xlims + ylims)

sp2.minorticks_on()
sp2.set_xlabel('(V - J)$_{\mathrm{rest}}$', size=22)
sp2.set_ylabel('(U - V)$_{\mathrm{rest}}$', size=22)
sp2.axis(xlims + ylims)


mypy.uvj_select_region(0.8, subplot=sp1, kolor='w', lw=8)
mypy.uvj_select_region(0.8, subplot=sp1, kolor='k', lw=3.5)
mypy.uvj_select_region(0.8, subplot=sp2, kolor='w', lw=8)
mypy.uvj_select_region(0.8, subplot=sp2, kolor='k', lw=3.5)







custom_colors = [pylab.cm.spectral(i) for i in pylab.linspace(0, 1, 20)]
custom_colors[0] = (1., 1., 1., 1.)
custom_colors[-1] = (0., 0., 0., 1.)

my_spectral = matplotlib.colors.ListedColormap(custom_colors, name='custom spectral_r', N=None)











vmin, vmax = 0, 0.5
implot = sp1.imshow(uv_errMap_zphot, interpolation='nearest', cmap=my_spectral,
	               vmin=vmin, vmax=vmax, extent=extent_zphot)

cbar = fig.colorbar(implot, ax=sp1, label='error  in  (U - V)$_{\mathrm{rest}}$')
cbar.set_clim(vmin, vmax)







vmin, vmax = 0, 0.75
implot = sp2.imshow(vj_errMap_zphot, interpolation='nearest', cmap=my_spectral,
	                vmin=vmin, vmax=vmax, extent=extent_zphot)

cbar = fig.colorbar(implot, ax=sp2, label='error  in  (V - J)$_{\mathrm{rest}}$')
cbar.set_clim(vmin, vmax)



pdfname = '../output/UVJ_Map_errors.pdf'
pylab.savefig(pdfname)
pylab.close()
print '\nwrote to:\n\n!open %s' % pdfname






















