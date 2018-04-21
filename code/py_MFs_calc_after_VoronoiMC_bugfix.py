
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
    if type(uv)==list or type(uv)==float or type(uv)==int or type(uv)==long or type(uv)==pylab.float64: uv = pylab.array(uv)
    if type(vj)==list or type(vj)==float or type(vj)==int or type(vj)==long or type(vj)==pylab.float64: vj = pylab.array(vj)

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


		###  UPDATING USE FLAG WITH REDUCE CHI**2 THRESHOLD
		chi2red = self.zout.chi_p / (self.zout.nfilt - 1.)
		cinds = pylab.find((chi2red > chi2red_thresh) & (self.cat.z_spec < 0))
		self.cat.use[cinds] = 0.


		###  UVJ classification
		self.UV_color = -2.5 * pylab.log10(self.restframe.restflux_U / self.restframe.restflux_V)
		self.VJ_color = -2.5 * pylab.log10(self.restframe.restflux_V / self.restframe.restflux_J)
		#self.uvj_class = mypy.uvj_select(self.UV_color, self.VJ_color, self.fout.z)
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
fields.append(field('SC1604',  'sc1604_v0.0.6',       0.910,  0.029, alpha=1./500, chi2red_thresh=10))
fields.append(field('SC0910',  'cl0910+5422_v0.0.5',  1.110,  0.035, alpha=1./500, chi2red_thresh=10))
fields.append(field('SC0849',  'sc0849+4452_v0.0.4',  1.261,  0.029, alpha=1./600, chi2red_thresh=10, uvj_intercept=0.54))


fields.append(field('CL1137',  'cl1137+3007_v0.0.2',  0.959,  0.032, alpha=1./600, chi2red_thresh=10))
fields.append(field('CL1350',  'cl1350+6007_v0.0.2',  0.804,  0.035, alpha=1./600, chi2red_thresh=10))
fields.append(field('CL1429',  'cl1429+4241_v0.0.4',  0.920,  0.051, alpha=1./500, chi2red_thresh=10))
fields.append(field('RXJ1053', 'rxj1053+5735_v0.0.3', 1.140,  0.031, alpha=1./600, chi2red_thresh=10))
fields.append(field('RXJ1221', 'rxj1221+4918_v0.0.4', 0.700,  0.023, alpha=1./600, chi2red_thresh=10))
fields.append(field('XLSS005', 'xlss005_v0.0.3',      1.000,  0.024, alpha=1./600, chi2red_thresh=10))
fields.append(field('SG0023',  'sg0023+0423_v0.2.0',  0.845,  0.025, alpha=1./600, chi2red_thresh=14))
print ''



fields[0].version_mugshots = 'nep200_v0.0.5'
fields[1].version_mugshots = 'sc1324_v0.0.2'
fields[2].version_mugshots = 'rcs0224-0002_v0.0.2'
fields[3].version_mugshots = 'rxj1716+6708_v0.0.9'
fields[4].version_mugshots = 'nep5281_v0.0.2'
fields[5].version_mugshots = 'sc1604_v0.0.6'
fields[6].version_mugshots = 'cl0910+5422_v0.0.3'
fields[7].version_mugshots = 'sc0849+4452_v0.0.2'
fields[8].version_mugshots = 'cl1137+3007_v0.0.1'
fields[9].version_mugshots = 'cl1350+6007_v0.0.1'
fields[10].version_mugshots = 'cl1429+4241_v0.0.2'
fields[11].version_mugshots = 'rxj1053+5735_v0.0.3'
fields[12].version_mugshots = 'rxj1221+4918_v0.0.2'
fields[13].version_mugshots = 'xlss005_v0.0.2'
fields[14].version_mugshots = 'sg0023+0423_v0.2.0'




for i in [12203, 13820]: fields[0].cat.use[i-1] *= 0
for i in [4386, 10388, 10907]: fields[1].cat.use[i-1] *= 0
for i in [17351, 26627, 38501, 43517, 51926, 55037, 55089, 61118]: fields[2].cat.use[i-1] *= 0
for i in []: fields[3].cat.use[i-1] *= 0
for i in [11331, 17375]: fields[4].cat.use[i-1] *= 0
for i in [105231, 170586]: fields[5].cat.use[i-1] *= 0
for i in [51114, 57254]: fields[6].cat.use[i-1] *= 0
for i in [9419, 15619, 18725, 22467, 25342, 28122]: fields[7].cat.use[i-1] *= 0



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







for f in fields:
	inds = pylab.find((f.cat.use == 1) &
		              (f.fout.z > 0.55) &
		              (f.fout.z < 1.30) &
		              (f.fout.lmass < 9.5+0.25/2) &
		              (f.fout.lmass > pylab.interp(f.fout.z, f.masslimits.z, f.masslimits.masslim_model_ssp)) &
		              (f.uvj_class == 0) & 
		              (f.voronoi.dens_at_z > 0.5))

	s = ''
	for id in f.cat.id[inds]:
		s += ' ./%s/mugshot%06i_%s.pdf' % (f.name, id, f.version_mugshots)

	print '\n%s\n' % s












t0 = time.time()

dm = 0.25
lmassbins = pylab.arange(9.5-dm/2., 11.75+dm, dm)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.
corr_factor = pylab.array([ 1.167, 1.091, 1.250, 0.917, 1.071, 1.231, 1.000, 1.000, 1.000, 1.000])
corr_factor_sf = pylab.array([ 0.895, 0.876, 0.787, 0.848, 0.756, 0.936, 0.791, 1.067, 1.000, 1.000])
corr_factor_qu = pylab.array([ 0.333, 1.455, 0.929, 0.831, 0.924, 0.892, 0.889, 0.858, 0.979, 0.979])
#corr_factor_sf = pylab.ones(len(corr_factor))
#corr_factor_qu = pylab.ones(len(corr_factor))



'''
##################################
###  calculating full MFs for LSS
##################################

schechter_fits = []
schechter_errs = []
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
	outer = open('../data_after_VoronoiMC_bugfix/MF_%s_LSS.dat' % f.name, 'w')
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
	print 'wrote to: ../data_after_VoronoiMC_bugfix/MF_%s_LSS.dat' % f.name






	###  UVJ class based on EAZY in 2-filter mode
	if True:

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
		outer = open('../data_after_VoronoiMC_bugfix/MF_%s_LSS_2filter.dat' % f.name, 'w')
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
		print 'wrote to: ../data_after_VoronoiMC_bugfix/MF_%s_LSS_2filter.dat' % f.name


print ''































######################################
###  calculating MFs in voronoi bins
######################################


###  correction factors for dm=0.25, lmass = pylab.arange(pylab.arange(9.5-dm/2., 11.5+dm, dm))
###    corr_factors[0]  -->  0.0 < log(1+delta) < 0.5
###    corr_factors[1]  -->  0.5 < log(1+delta) < 1.0
###    corr_factors[2]  -->  1.0 < log(1+delta) < 1.5
###    corr_factors[3]  -->  1.5 < log(1+delta) < 2.0

dm = 0.25
lmassbins = pylab.arange(9.5-dm/2., 11.75+dm, dm)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.

corr_factors = pylab.array([[ 0.786, 0.738, 0.676, 0.713, 0.632, 0.769, 0.616, 0.744, 0.786, 0.786], 
	                        [ 0.957, 1.077, 0.940, 0.918, 0.906, 0.974, 0.941, 0.978, 1.125, 1.125], 
	                        [ 1.000, 1.294, 1.000, 1.029, 1.066, 1.020, 1.105, 0.906, 1.062, 1.062], 
	                        [ 1.167, 1.091, 1.250, 0.917, 1.071, 1.231, 1.000, 1.000, 1.000, 1.000]])

elo_corrfactos = pylab.array([[ 0.103, 0.095, 0.104, 0.081, 0.078, 0.059, 0.087, 0.118, 0.149, 0.149], 
	                          [ 0.164, 0.186, 0.103, 0.081, 0.072, 0.064, 0.057, 0.064, 0.191, 0.191], 
	                          [ 0.300, 0.189, 0.134, 0.092, 0.060, 0.074, 0.061, 0.065, 0.052, 0.052], 
	                          [ 0.138, 0.139, 0.228, 0.054, 0.082, 0.091, 0.033, 0.000, 0.000, 0.000]])

ehi_corrfactos = pylab.array([[ 0.131, 0.119, 0.125, 0.099, 0.093, 0.072, 0.105, 0.159, 0.278, 0.278], 
	                          [ 0.237, 0.248, 0.135, 0.107, 0.097, 0.085, 0.082, 0.103, 0.435, 0.435], 
	                          [ 0.496, 0.288, 0.202, 0.140, 0.087, 0.107, 0.107, 0.122, 0.183, 0.183], 
	                          [ 0.488, 0.317, 0.463, 0.133, 0.147, 0.154, 0.092, 0.172, 0.199, 0.199]])



corr_factors_sf = pylab.array([[ 0.815, 0.698, 0.673, 0.679, 0.549, 0.788, 0.571, 1.000, 0.667, 0.667], 
	                           [ 0.955, 1.040, 0.892, 1.000, 0.848, 0.949, 0.906, 1.154, 1.500, 1.500], 
	                           [ 1.000, 1.333, 0.944, 1.067, 1.067, 1.136, 1.111, 0.875, 1.000, 1.000], 
	                           [ 1.167, 1.000, 1.250, 0.857, 1.200, 1.667, 1.000, 1.000, 1.000, 1.000]])

elo_corrfactos_sf = pylab.array([[ 0.104, 0.096, 0.122, 0.109, 0.115, 0.099, 0.141, 0.329, 0.276, 0.276], 
	                             [ 0.172, 0.189, 0.130, 0.105, 0.114, 0.107, 0.078, 0.140, 0.414, 0.414], 
	                             [ 0.300, 0.225, 0.178, 0.214, 0.077, 0.159, 0.170, 0.104, 0.000, 0.000], 
	                             [ 0.138, 0.146, 0.383, 0.092, 0.306, 0.318, 0.167, 0.000, 0.000, 0.000]])

ehi_corrfactos_sf = pylab.array([[ 0.133, 0.120, 0.150, 0.137, 0.144, 0.129, 0.185, 0.588, 0.976, 0.976], 
	                             [ 0.248, 0.254, 0.174, 0.143, 0.165, 0.149, 0.134, 0.285, 1.463, 1.463], 
	                             [ 0.496, 0.370, 0.273, 0.328, 0.137, 0.236, 0.387, 0.366, 0.646, 0.646], 
	                             [ 0.488, 0.404, 0.871, 0.229, 0.697, 0.607, 0.462, 0.000, 2.582, 2.582]])



corr_factors_qu = pylab.array([[ 1.000, 2.000, 0.688, 0.774, 0.727, 0.754, 0.647, 0.688, 0.818, 0.818], 
	                           [ 1.000, 2.000, 1.077, 0.737, 0.968, 1.000, 0.972, 0.906, 1.000, 1.000], 
	                           [ 1.000, 1.200, 1.143, 1.000, 1.065, 0.931, 1.103, 0.917, 1.083, 1.083], 
	                           [ 1.000, 1.333, 1.250, 1.000, 1.043, 1.100, 1.000, 1.000, 1.000, 1.000]])

elo_corrfactos_qu = pylab.array([[ 0.644, 0.644, 0.180, 0.113, 0.100, 0.071, 0.110, 0.122, 0.166, 0.166], 
	                             [ 0.000, 0.829, 0.118, 0.113, 0.081, 0.069, 0.080, 0.065, 0.195, 0.195], 
	                             [ 0.000, 0.306, 0.118, 0.000, 0.087, 0.044, 0.056, 0.076, 0.069, 0.069], 
	                             [ 0.000, 0.276, 0.207, 0.000, 0.067, 0.064, 0.000, 0.000, 0.000, 0.000]])

ehi_corrfactos_qu = pylab.array([[ 1.600, 1.600, 0.287, 0.166, 0.136, 0.094, 0.140, 0.173, 0.336, 0.336], 
	                             [ 2.582, 2.927, 0.268, 0.202, 0.138, 0.114, 0.128, 0.122, 0.539, 0.539], 
	                             [ 0.000, 0.697, 0.418, 0.129, 0.144, 0.110, 0.118, 0.154, 0.244, 0.244], 
	                             [ 0.000, 0.976, 0.732, 0.258, 0.151, 0.160, 0.092, 0.172, 0.215, 0.215]])


for i in range(len(fields)):

	f = fields[i]


	###  finding lower/upper zphot limits
	zlo_phot = f.zclust - 1.5 * f.sigmaz * (1 + f.zclust)
	zhi_phot = f.zclust + 1.5 * f.sigmaz * (1 + f.zclust)


	###  finding lower/upper zspec limits
	zlo_spec, zhi_spec = f.zclust, f.zclust
	for ss in f.substructures:
		vmean = vrecess(ss.zmean)
		vlo = vmean - 3*ss.vsigma1mpc
		vhi = vmean + 3*ss.vsigma1mpc
		zlo = zrecess(vlo)
		zhi = zrecess(vhi)
		if zlo < zlo_spec: zlo_spec = zlo
		if zhi > zhi_spec: zhi_spec = zhi









	#print '\nUSING EXPANDED ZSPEC WINDOW !!!\n'
	#zlo_spec = f.zclust - 0.2
	#zhi_spec = f.zclust + 0.2







	###  LSS members selected by zspec
	members_zspec = pylab.find((f.cat.use[f.inds_spatial] == 1) &
		                       (-pylab.isnan(f.fout.lmass[f.inds_spatial])) &
		                       (f.cat.z_spec[f.inds_spatial] > zlo_spec) &
    		                   (f.cat.z_spec[f.inds_spatial] < zhi_spec))

	###  LSS candidates selected by zphot
	members_zphot = pylab.find((f.cat.use[f.inds_spatial] == 1) &
		                       (-pylab.isnan(f.fout.lmass[f.inds_spatial])) &
		                       (f.cat.z_spec[f.inds_spatial] < 0) &
		                       (f.fout.z[f.inds_spatial] > zlo_phot) &
   			                   (f.fout.z[f.inds_spatial] < zhi_phot))



	###  assigning Voroini overdensity values to galaxies
	overdens_zspec = []
	overdens_zphot = []

	###  ... for zspec members this is easy, grab the zslice that the zspec is inside
	for member in members_zspec:
		zdif = abs(f.cat.z_spec[f.inds_spatial][member] - f.voronoi.zbars)
		vi = pylab.find(zdif == zdif.min())[0]
		overdens_zspec.append(f.voronoi.overdens_matrix[f.inds_spatial][member][vi])

	###  ... for zphot it's a bit trickier. What I will do is adopt the zmean of
	###      the nearest substructure in projection.
	for member in members_zphot:
		ra0, dec0 = f.cat.ra[f.inds_spatial][member], f.cat.dec[f.inds_spatial][member]
		rproj_list = []
		for ss in f.substructures:
			rproj_arcsec = mypy.radec_sep(ra0, dec0, ss.ra, ss.dec)
			rproj_list.append(rproj_arcsec)

		rind = rproj_list.index(min(rproj_list))
		zdif = abs(f.substructures[rind].zmean - f.voronoi.zbars)
		vi = pylab.find(zdif == zdif.min())[0]
		overdens_zphot.append(f.voronoi.overdens_matrix[f.inds_spatial][member][vi])

	overdens_zspec = pylab.array(overdens_zspec)
	overdens_zphot = pylab.array(overdens_zphot)




	overdens_bins = pylab.array([0., 0.5, 1., 1.5, 2.])
	overdens_bars = (overdens_bins[1:] + overdens_bins[:-1]) / 2.
	overdens_labels = ['00v05', '05v10', '10v15', '15v20']

	###  binning galaxies by overdensity
	digi_overdens_zspec = pylab.digitize(overdens_zspec, overdens_bins)
	digi_overdens_zphot = pylab.digitize(overdens_zphot, overdens_bins)

	###  binning galaxies by stellar mass
	digi_lmass_zspec = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zspec], lmassbins)
	digi_lmass_zphot = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zphot], lmassbins)



	###  grabbing UVJ class
	inds_sf_zspec = pylab.find(f.uvj_class[f.inds_spatial][members_zspec] == 1)
	inds_sf_zphot = pylab.find(f.uvj_class[f.inds_spatial][members_zphot] == 1)
	inds_qu_zspec = pylab.find(f.uvj_class[f.inds_spatial][members_zspec] == 0)
	inds_qu_zphot = pylab.find(f.uvj_class[f.inds_spatial][members_zphot] == 0)

	digi_overdens_zspec_sf = pylab.digitize(overdens_zspec[inds_sf_zspec], overdens_bins)
	digi_overdens_zphot_sf = pylab.digitize(overdens_zphot[inds_sf_zphot], overdens_bins)
	digi_lmass_zspec_sf = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zspec][inds_sf_zspec], lmassbins)
	digi_lmass_zphot_sf = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zphot][inds_sf_zphot], lmassbins)

	digi_overdens_zspec_qu = pylab.digitize(overdens_zspec[inds_qu_zspec], overdens_bins)
	digi_overdens_zphot_qu = pylab.digitize(overdens_zphot[inds_qu_zphot], overdens_bins)
	digi_lmass_zspec_qu = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zspec][inds_qu_zspec], lmassbins)
	digi_lmass_zphot_qu = pylab.digitize(f.fout.lmass[f.inds_spatial][members_zphot][inds_qu_zphot], lmassbins)




	for overi in range(1, len(overdens_bins)):

		###  writing output
		outer = open('../data_after_VoronoiMC_bugfix/MF_%s_LSS_%s.dat' % (f.name, overdens_labels[overi-1]), 'w')
		outer.write('#   Ntot = Nzspec + Nzphot_corr\n')
		outer.write('#   Nlo_tot = lower CI for Ntot (nlo_spec**2 + nlo_phot**2 + corr_f**2)**0.5\n')
		outer.write('#   Nhi_tot = upper CI for Ntot (nhi_spec**2 + nhi_phot**2 + corr_f**2)**0.5\n')
		outer.write('#   Nzspec = number of LSS members selected by zspec\n')
		outer.write('#   Nzphot = number of candidate LSS members selected in broad zphot bin\n')
		outer.write('#   Nzphot_corr = Nzphot corrected by false +/- fractions to estimate the\n')
		outer.write('#                 number of LSS members missed by spectroscopy\n')
		outer.write('# lmass Ntot Nlo_tot Nhi_tot Nzspec Nzphot Nzphot_corr')
		outer.write('   Nsf Nzspec_sf Nlo_sf Nhi_sf Nzphot_sf Nzphot_corr_sf')
		outer.write('   Nqu Nzspec_qu Nlo_qu Nhi_qu Nzphot_qu Nzphot_corr_qu\n')


		for mi in range(1, len(lmassbins)):


			###################
			###  ALL GALAXIES
			###################
			count_zspec = len(pylab.find((digi_overdens_zspec == overi) & (digi_lmass_zspec == mi)))
			count_zphot = len(pylab.find((digi_overdens_zphot == overi) & (digi_lmass_zphot == mi)))
			count_zphot_corr = count_zphot * corr_factors[overi-1][mi-1]
			nfinal = count_zspec + count_zphot_corr


			###  propigating errors
			count_hi_zspec, count_lo_zspec = mypy.massfunc.confidence_interval(count_zspec)
			count_hi_zphot_corr, count_lo_zphot_corr = mypy.massfunc.confidence_interval(count_zphot_corr)

			frac_elo = elo_corrfactos[overi-1][mi-1] / corr_factors[overi-1][mi-1]
			frac_ehi = ehi_corrfactos[overi-1][mi-1] / corr_factors[overi-1][mi-1]

			count_lo_zphot_corr = pylab.sqrt(count_lo_zphot_corr**2 + count_zphot_corr**2 * frac_elo**2)
			count_hi_zphot_corr = pylab.sqrt(count_hi_zphot_corr**2 + count_zphot_corr**2 * frac_ehi**2)

			elo = pylab.sqrt(count_lo_zspec**2 + count_lo_zphot_corr**2)
			ehi = pylab.sqrt(count_hi_zspec**2 + count_hi_zphot_corr**2)


			outer.write(' %5.2f' % lmassbars[mi-1])
			outer.write(' %4i' % (nfinal + 0.5))
			outer.write(' %4i' % elo)
			outer.write(' %4i' % ehi)
			outer.write(' %4i' % count_zspec)
			outer.write(' %4i' % count_zphot)
			outer.write(' %4i' % count_zphot_corr)




			###################
			###  SF GALAXIES
			###################
			count_zspec_sf = len(pylab.find((digi_overdens_zspec_sf == overi) & (digi_lmass_zspec_sf == mi)))
			count_zphot_sf = len(pylab.find((digi_overdens_zphot_sf == overi) & (digi_lmass_zphot_sf == mi)))
			count_zphot_corr_sf = count_zphot_sf * corr_factors_sf[overi-1][mi-1]
			nfinal_sf = count_zspec_sf + count_zphot_corr_sf


			###  propigating errors
			count_hi_zspec_sf, count_lo_zspec_sf = mypy.massfunc.confidence_interval(count_zspec_sf)
			count_hi_zphot_corr_sf, count_lo_zphot_corr_sf = mypy.massfunc.confidence_interval(count_zphot_corr_sf)

			frac_elo_sf = elo_corrfactos_sf[overi-1][mi-1] / corr_factors_sf[overi-1][mi-1]
			frac_ehi_sf = ehi_corrfactos_sf[overi-1][mi-1] / corr_factors_sf[overi-1][mi-1]

			count_lo_zphot_corr_sf = pylab.sqrt(count_lo_zphot_corr_sf**2 + count_zphot_corr_sf**2 * frac_elo_sf**2)
			count_hi_zphot_corr_sf = pylab.sqrt(count_hi_zphot_corr_sf**2 + count_zphot_corr_sf**2 * frac_ehi_sf**2)

			elo_sf = pylab.sqrt(count_lo_zspec_sf**2 + count_lo_zphot_corr_sf**2)
			ehi_sf = pylab.sqrt(count_hi_zspec_sf**2 + count_hi_zphot_corr_sf**2)


			outer.write('     ')
			outer.write(' %4i' % (nfinal_sf + 0.5))
			outer.write(' %4i' % elo_sf)
			outer.write(' %4i' % ehi_sf)
			outer.write(' %4i' % count_zspec_sf)
			outer.write(' %4i' % count_zphot_sf)
			outer.write(' %4i' % count_zphot_corr_sf)




			###################
			###  QU GALAXIES
			###################
			count_zspec_qu = len(pylab.find((digi_overdens_zspec_qu == overi) & (digi_lmass_zspec_qu == mi)))
			count_zphot_qu = len(pylab.find((digi_overdens_zphot_qu == overi) & (digi_lmass_zphot_qu == mi)))
			count_zphot_corr_qu = count_zphot_qu * corr_factors_qu[overi-1][mi-1]
			nfinal_qu = count_zspec_qu + count_zphot_corr_qu


			###  propigating errors
			count_hi_zspec_qu, count_lo_zspec_qu = mypy.massfunc.confidence_interval(count_zspec_qu)
			count_hi_zphot_corr_qu, count_lo_zphot_corr_qu = mypy.massfunc.confidence_interval(count_zphot_corr_qu)

			frac_elo_qu = elo_corrfactos_qu[overi-1][mi-1] / corr_factors_qu[overi-1][mi-1]
			frac_ehi_qu = ehi_corrfactos_qu[overi-1][mi-1] / corr_factors_qu[overi-1][mi-1]

			count_lo_zphot_corr_qu = pylab.sqrt(count_lo_zphot_corr_qu**2 + count_zphot_corr_qu**2 * frac_elo_qu**2)
			count_hi_zphot_corr_qu = pylab.sqrt(count_hi_zphot_corr_qu**2 + count_zphot_corr_qu**2 * frac_ehi_qu**2)

			elo_qu = pylab.sqrt(count_lo_zspec_qu**2 + count_lo_zphot_corr_qu**2)
			ehi_qu = pylab.sqrt(count_hi_zspec_qu**2 + count_hi_zphot_corr_qu**2)


			outer.write('     ')
			outer.write(' %4i' % (nfinal_qu + 0.5))
			outer.write(' %4i' % elo_qu)
			outer.write(' %4i' % ehi_qu)
			outer.write(' %4i' % count_zspec_qu)
			outer.write(' %4i' % count_zphot_qu)
			outer.write(' %4i' % count_zphot_corr_qu)


			outer.write('\n')
		outer.close()
		print 'wrote to: ../data_after_VoronoiMC_bugfix/MF_%s_LSS_%s.dat' % (f.name, overdens_labels[overi-1])
	print ''












#####################################
###  calculating full MFs for field
#####################################

class massFunction_field:
	def __init__(self, name, area, lmassbars, zlos, zhis):
		self.name = name
		self.area = area    # unobscured survey area in arcmin**2
		self.lmassbars = lmassbars
		self.zlos = zlos
		self.zhis = zhis
		self.zbars = (zhis + zlos) / 2.
		self.dz = (zhis - zlos)
		self.number_counts = pylab.array([pylab.zeros(len(lmassbars)) for zi in range(len(zlos))])
		self.nlo_poisson = pylab.array([pylab.zeros(len(lmassbars)) for zi in range(len(zlos))])
		self.nhi_poisson = pylab.array([pylab.zeros(len(lmassbars)) for zi in range(len(zlos))])

dm = 0.25
lmassbins = pylab.arange(8.5-dm/2., 11.75+dm, dm)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2.

zbins = pylab.array([[0.20, 0.50],
                     [0.50, 0.75],
                     [0.75, 1.00],
                     [1.00, 1.25],
                     [1.25, 1.50],
                     [1.50, 2.00]])



outer_massiveGals = open('../output_after_VoronoiMC_bugfix/massiveGals.txt' ,'w')
for fi in range(len(fields)):

	f = fields[fi]
	mf_field = massFunction_field(f.name, f.area_arcmin2, lmassbars, zbins[:,0], zbins[:,1])
	idx = []

	for zi in range(len(zbins)):

		zlo, zhi = zbins[zi]

		subinds = pylab.find((f.cat.use[f.inds_spatial] == 1) &
			                 (f.zout.odds[f.inds_spatial] > 0.75) &
			                 (f.fout.z[f.inds_spatial] > zlo) &
			                 (f.fout.z[f.inds_spatial] < zhi))

		subinds_massive = pylab.find((f.cat.use[f.inds_spatial] == 1) &
					                 (f.zout.odds[f.inds_spatial] > 0.75) &
			                         (f.fout.z[f.inds_spatial] > zlo) &
			                         (f.fout.z[f.inds_spatial] < zhi) &
			                         (f.fout.lmass[f.inds_spatial] >= 11.25))

		idx += pylab.sort(f.cat.id[f.inds_spatial][subinds_massive]).tolist()


		digi_mass = pylab.digitize(f.fout.lmass[f.inds_spatial][subinds], lmassbins)

		ngal_bins = pylab.bincount(digi_mass, minlength=len(lmassbins)+1)[1:-1]
		nlo_poisson, nhi_poisson = [], []
		for n in ngal_bins:
			nhi, nlo = mypy.massfunc.confidence_interval(n)
			nlo_poisson.append(nlo)
			nhi_poisson.append(nhi)
		nlo_poisson, nhi_poisson = pylab.array(nlo_poisson), pylab.array(nhi_poisson)


		mf_field.number_counts[zi] += ngal_bins
		mf_field.nlo_poisson[zi] += nlo_poisson
		mf_field.nhi_poisson[zi] += nhi_poisson

	pickle.dump(mf_field, open('../data_after_VoronoiMC_bugfix/MF_%s_field.pickle' % f.name, 'wb'))
	print 'wrote to: ../data_after_VoronoiMC_bugfix/MF_%s_field.pickle' % f.name


	for idx_i in idx:
		outer_massiveGals.write(' %s/mugshot%05i_%s.pdf' % (f.version, idx_i, f.version))

outer_massiveGals.close()
print ''


tf = time.time()
print '%.1f mins elapsed\n\n' % ((tf-t0)/60.)
'''













################################################
###  calculating MFs for each Voronoi zslice
###  NOTE: these MFs do not have a statistical
###     correction applied. Just number counts
################################################

from scipy import interpolate
#voronoi_dir = '/Users/atomczak/DATA/ORELSE/Voronoi'
voronoi_dir = '/Users/atomczak/GoogleDrive/ORELSE/Voronoi_Maps/Number_Density_and_Overdensity_Maps'

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
		self.lmassbins = pylab.arange(8.0-self.dm/2., 11.75+self.dm, self.dm)
		self.lmassbars = (self.lmassbins[1:] + self.lmassbins[:-1]) / 2.

		self.areas_full_arcmin2 = pylab.zeros(len(self.zbars))               # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_field_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_v_lt_0_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_n10v05_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_n05v00_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_00v05_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_05v10_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_10v15_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.areas_15v20_arcmin2 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice

		self.volumes_full_Mpc3 = pylab.zeros(len(self.zbars))               # volumes corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_field_Mpc3 = pylab.zeros(len(self.zbars))              # volumes corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_v_lt_0_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_n10v05_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
		self.volumes_n05v00_Mpc3 = pylab.zeros(len(self.zbars))              # areas corresponding to log(overdensity)<0.3 for each zslice
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
		self.overdens_bins = pylab.array([-1., -0.5, 0., 0.5, 1., 1.5, 2.])
		self.overdens_bars = (self.overdens_bins[1:] + self.overdens_bins[:-1]) / 2.

		self.ngals_v_lt_0 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_v_lt_0 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_v_lt_0 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_v_lt_0 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n10v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n10v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n10v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n10v05 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n05v00 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n05v00 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n05v00 = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n05v00 = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

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

		self.areas_all_arcmin2 = [self.areas_n10v05_arcmin2, self.areas_n05v00_arcmin2, self.areas_00v05_arcmin2, self.areas_05v10_arcmin2, self.areas_10v15_arcmin2, self.areas_15v20_arcmin2]
		self.volumes_all_Mpc3 = [self.volumes_n10v05_Mpc3, self.volumes_n05v00_Mpc3, self.volumes_00v05_Mpc3, self.volumes_05v10_Mpc3, self.volumes_10v15_Mpc3, self.volumes_15v20_Mpc3]
		self.ngals_all = [self.ngals_n10v05, self.ngals_n05v00, self.ngals_00v05, self.ngals_05v10, self.ngals_10v15, self.ngals_15v20]
		self.elo_ngals_all = [self.elo_ngals_n10v05, self.elo_ngals_n05v00, self.elo_ngals_00v05, self.elo_ngals_05v10, self.elo_ngals_10v15, self.elo_ngals_15v20]
		self.ehi_ngals_all = [self.ehi_ngals_n10v05, self.ehi_ngals_n05v00, self.ehi_ngals_00v05, self.ehi_ngals_05v10, self.ehi_ngals_10v15, self.ehi_ngals_15v20]
		self.phis_all = [self.phis_n10v05, self.phis_n05v00, self.phis_00v05, self.phis_05v10, self.phis_10v15, self.phis_15v20]



		##################
		###  SF galaxies
		##################

		self.ngals_v_lt_0_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_v_lt_0_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_v_lt_0_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_v_lt_0_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n10v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n10v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n10v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n10v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n05v00_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n05v00_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n05v00_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n05v00_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_full_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_full_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_full_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_full_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_field_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_field_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_field_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_field_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_00v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_00v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_00v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_00v05_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_05v10_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_05v10_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_05v10_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_05v10_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_10v15_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_10v15_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_10v15_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_10v15_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_15v20_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_15v20_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_15v20_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_15v20_sf = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_all_sf = [self.ngals_n10v05_sf, self.ngals_n05v00_sf, self.ngals_00v05_sf, self.ngals_05v10_sf, self.ngals_10v15_sf, self.ngals_15v20_sf]
		self.elo_ngals_all_sf = [self.elo_ngals_n10v05_sf, self.elo_ngals_n05v00_sf, self.elo_ngals_00v05_sf, self.elo_ngals_05v10_sf, self.elo_ngals_10v15_sf, self.elo_ngals_15v20_sf]
		self.ehi_ngals_all_sf = [self.ehi_ngals_n10v05_sf, self.ehi_ngals_n05v00_sf, self.ehi_ngals_00v05_sf, self.ehi_ngals_05v10_sf, self.ehi_ngals_10v15_sf, self.ehi_ngals_15v20_sf]
		self.phis_all_sf = [self.phis_n10v05_sf, self.phis_n05v00_sf, self.phis_00v05_sf, self.phis_05v10_sf, self.phis_10v15_sf, self.phis_15v20_sf]



		##################
		###  QU galaxies
		##################

		self.ngals_full_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_full_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_full_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_full_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_field_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_field_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_field_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_field_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_v_lt_0_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_v_lt_0_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_v_lt_0_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_v_lt_0_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n10v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n10v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n10v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n10v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_n05v00_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_n05v00_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_n05v00_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_n05v00_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_00v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_00v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_00v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_00v05_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_05v10_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_05v10_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_05v10_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_05v10_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_10v15_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_10v15_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_10v15_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_10v15_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_15v20_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # number of galaxies per massbin for each zslice
		self.elo_ngals_15v20_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # lower CI for Poisson stats
		self.ehi_ngals_15v20_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))  # upper CI for Poisson stats
		self.phis_15v20_qu = pylab.zeros((len(self.zbars), len(self.lmassbars)))      # volume density of galaxies per massbin for each zslice

		self.ngals_all_qu = [self.ngals_n10v05_qu, self.ngals_n05v00_qu, self.ngals_00v05_qu, self.ngals_05v10_qu, self.ngals_10v15_qu, self.ngals_15v20_qu]
		self.elo_ngals_all_qu = [self.elo_ngals_n10v05_qu, self.elo_ngals_n05v00_qu, self.elo_ngals_00v05_qu, self.elo_ngals_05v10_qu, self.elo_ngals_10v15_qu, self.elo_ngals_15v20_qu]
		self.ehi_ngals_all_qu = [self.ehi_ngals_n10v05_qu, self.ehi_ngals_n05v00_qu, self.ehi_ngals_00v05_qu, self.ehi_ngals_05v10_qu, self.ehi_ngals_10v15_qu, self.ehi_ngals_15v20_qu]
		self.phis_all_qu = [self.phis_n10v05_qu, self.phis_n05v00_qu, self.phis_00v05_qu, self.phis_05v10_qu, self.phis_10v15_qu, self.phis_15v20_qu]



images = ['/Users/atomczak/GoogleDrive/ORELSE/images/nep200_v0.0.6/nep200_v0.0.6_detection_ri.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/sc1324_v0.0.3/sc1324_v0.0.3_detection_i.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/rcs0224-0002_v0.0.3/rcs0224-0002_v0.0.3_detection_I+.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/rxj1716+6708_v0.0.8/rxj1716+6708_v0.0.8_detection_RcI+Z+.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/nep5281_v0.0.3/nep5281_v0.0.3_detection_Y.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/sc1604_v0.0.4/sc1604_v0.0.4_detection_Rc.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/cl0910+5422_v0.0.4/cl0910+5422_v0.0.4_detection_RcI+Z+.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/sc0849+4452_v0.0.3/sc0849+4452_v0.0.3_detection_Z+.fits']

images = ['/Users/atomczak/GoogleDrive/ORELSE/images/cl1137+3007_v0.0.1/cl1137+3007_v0.0.1_detection_Z+.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/cl1350+6007_v0.0.2/cl1350+6007_v0.0.2_detection_r_mega.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/cl1429+4241_v0.0.3/cl1429+4241_v0.0.3_detection_Y.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/rxj1053+5735_v0.0.2/rxj1053+5735_v0.0.2_detection_Z+.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/rxj1221+4918_v0.0.2/rxj1221+4918_v0.0.2_detection_i.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/xlss005_v0.0.3/xlss005_v0.0.3_detection_i_mega.fits',
'/Users/atomczak/GoogleDrive/ORELSE/images/sg0023+0423_v0.2.0/sg0023+0423_v0.2.0_detection_R+I+.fits']




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
		inds_v_lt_0 = pylab.where((segmap_vor == 1) & (vmap <= 0))
		mf_voronoi_slices_field.areas_full_arcmin2[j] = len(inds_full[0]) * (px_scale_vor/60.)**2
		mf_voronoi_slices_field.areas_field_arcmin2[j] = len(inds_field[0]) * (px_scale_vor/60.)**2
		mf_voronoi_slices_field.areas_v_lt_0_arcmin2[j] = len(inds_v_lt_0[0]) * (px_scale_vor/60.)**2
		mf_voronoi_slices_field.volumes_full_Mpc3[j] = mypy.vcomoving_slice(mf_voronoi_slices_field.areas_full_arcmin2[j], zlo, zhi)
		mf_voronoi_slices_field.volumes_field_Mpc3[j] = mypy.vcomoving_slice(mf_voronoi_slices_field.areas_field_arcmin2[j], zlo, zhi)
		mf_voronoi_slices_field.volumes_v_lt_0_Mpc3[j] = mypy.vcomoving_slice(mf_voronoi_slices_field.areas_v_lt_0_arcmin2[j], zlo, zhi)

		vmap_vbin = vmap * 0.
		vmap_vbin[inds_full] = vmap[inds_full]
		#fits.writeto('../data_after_VoronoiMC_bugfix/voronoi_maps/%s_full_z%.4f.fits' % (f.name, mf_voronoi_slices_field.zbars[j]),
		#             vmap_vbin, header=vheader, clobber=1)
		vmap_vbin = vmap * 0.
		vmap_vbin[inds_field] = vmap[inds_field]
		#fits.writeto('../data_after_VoronoiMC_bugfix/voronoi_maps/%s_field_z%.4f.fits' % (f.name, mf_voronoi_slices_field.zbars[j]),
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
			#fits.writeto('../data_after_VoronoiMC_bugfix/voronoi_maps/%s_%.1fv%.1f_z%.4f.fits' % (f.name, mf_voronoi_slices_field.overdens_bins[k], mf_voronoi_slices_field.overdens_bins[k+1], mf_voronoi_slices_field.zbars[j]),
			#             vmap_vbin, header=vheader, clobber=1)





		###  number counting galaxies in massbins
		inds_full = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                   (f.fout.z[f.inds_spatial] >= zlo) &
			                   (f.fout.z[f.inds_spatial] < zhi))
		inds_field = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                    (f.fout.z[f.inds_spatial] >= zlo) &
			                    (f.fout.z[f.inds_spatial] < zhi) &
			                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= 0.3))
		inds_v_lt_0 = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                    (f.fout.z[f.inds_spatial] >= zlo) &
			                    (f.fout.z[f.inds_spatial] < zhi) &
			                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= 0))

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

		if len(inds_v_lt_0) > 0:
			digi_mass_v_lt_0 = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_v_lt_0], mf_voronoi_slices_field.lmassbins)
			ngal_v_lt_0 = pylab.bincount(digi_mass_v_lt_0, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
			nlo_v_lt_0, nhi_v_lt_0 = mypy.massfunc.CI(ngal_v_lt_0)
			mf_voronoi_slices_field.ngals_v_lt_0[j] += ngal_v_lt_0
			mf_voronoi_slices_field.elo_ngals_v_lt_0[j] += nlo_v_lt_0
			mf_voronoi_slices_field.ehi_ngals_v_lt_0[j] += nhi_v_lt_0
			mf_voronoi_slices_field.phis_v_lt_0[j] += ngal_v_lt_0 / mf_voronoi_slices_field.volumes_v_lt_0_Mpc3[j] / mf_voronoi_slices_field.dm
		else:
			nlo_v_lt_0, nhi_v_lt_0 = mypy.massfunc.CI(ngal_v_lt_0)
			mf_voronoi_slices_field.elo_ngals_v_lt_0[j] += nlo_v_lt_0
			mf_voronoi_slices_field.ehi_ngals_v_lt_0[j] += nhi_v_lt_0


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




		###############################
		###  repeating for SF galaxies
		###############################

		###  number counting galaxies in massbins
		inds_full = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                   (f.fout.z[f.inds_spatial] >= zlo) &
			                   (f.fout.z[f.inds_spatial] < zhi) &
			                   (f.uvj_class[f.inds_spatial] == 1))
		inds_field = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                    (f.fout.z[f.inds_spatial] >= zlo) &
			                    (f.fout.z[f.inds_spatial] < zhi) &
			                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= 0.3) &
			                    (f.uvj_class[f.inds_spatial] == 1))
		inds_v_lt_0 = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                    (f.fout.z[f.inds_spatial] >= zlo) &
			                    (f.fout.z[f.inds_spatial] < zhi) &
			                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= 0) &
			                    (f.uvj_class[f.inds_spatial] == 1))

		digi_mass_full = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_full], mf_voronoi_slices_field.lmassbins)
		ngal_full = pylab.bincount(digi_mass_full, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
		nlo_full, nhi_full = mypy.massfunc.CI(ngal_full)
		mf_voronoi_slices_field.ngals_full_sf[j] += ngal_full
		mf_voronoi_slices_field.elo_ngals_full_sf[j] += nlo_full
		mf_voronoi_slices_field.ehi_ngals_full_sf[j] += nhi_full
		mf_voronoi_slices_field.phis_full_sf[j] += ngal_full / mf_voronoi_slices_field.volumes_full_Mpc3[j] / mf_voronoi_slices_field.dm

		digi_mass_field = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_field], mf_voronoi_slices_field.lmassbins)
		ngal_field = pylab.bincount(digi_mass_field, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
		nlo_field, nhi_field = mypy.massfunc.CI(ngal_field)
		mf_voronoi_slices_field.ngals_field_sf[j] += ngal_field
		mf_voronoi_slices_field.elo_ngals_field_sf[j] += nlo_field
		mf_voronoi_slices_field.ehi_ngals_field_sf[j] += nhi_field
		mf_voronoi_slices_field.phis_field_sf[j] += ngal_field / mf_voronoi_slices_field.volumes_field_Mpc3[j] / mf_voronoi_slices_field.dm

		if len(inds_v_lt_0) > 0:
			digi_mass_v_lt_0 = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_v_lt_0], mf_voronoi_slices_field.lmassbins)
			ngal_v_lt_0 = pylab.bincount(digi_mass_v_lt_0, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
			nlo_v_lt_0, nhi_v_lt_0 = mypy.massfunc.CI(ngal_v_lt_0)
			mf_voronoi_slices_field.ngals_v_lt_0_sf[j] += ngal_v_lt_0
			mf_voronoi_slices_field.elo_ngals_v_lt_0_sf[j] += nlo_v_lt_0
			mf_voronoi_slices_field.ehi_ngals_v_lt_0_sf[j] += nhi_v_lt_0
			mf_voronoi_slices_field.phis_v_lt_0_sf[j] += ngal_v_lt_0 / mf_voronoi_slices_field.volumes_v_lt_0_Mpc3[j] / mf_voronoi_slices_field.dm
		else:
			nlo_v_lt_0, nhi_v_lt_0 = mypy.massfunc.CI(ngal_v_lt_0)
			mf_voronoi_slices_field.elo_ngals_v_lt_0_sf[j] += nlo_v_lt_0
			mf_voronoi_slices_field.ehi_ngals_v_lt_0_sf[j] += nhi_v_lt_0

		###  ... for overdensity bins
		for k in range(len(mf_voronoi_slices_field.overdens_bars)):
			inds_vbin =  pylab.find((f.cat.use[f.inds_spatial] == 1) & 
				                    (f.fout.z[f.inds_spatial] >= zlo) &
				                    (f.fout.z[f.inds_spatial] < zhi) &
				                    (f.uvj_class[f.inds_spatial] == 1) &
				                    (f.voronoi.overdens_matrix[f.inds_spatial, j] >= mf_voronoi_slices_field.overdens_bins[k]) &
				                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= mf_voronoi_slices_field.overdens_bins[k+1]))

			###  if no galaxies exist in this overdensity bin do nothing
			if len(inds_vbin) == 0:
				pass
			else:
				digi_mass_vbin = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_vbin], mf_voronoi_slices_field.lmassbins)
				ngal_vbin = pylab.bincount(digi_mass_vbin, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
				nlo_vbin, nhi_vbin = mypy.massfunc.CI(ngal_vbin)
				mf_voronoi_slices_field.ngals_all_sf[k][j] += ngal_vbin
				mf_voronoi_slices_field.elo_ngals_all_sf[k][j] += nlo_vbin
				mf_voronoi_slices_field.ehi_ngals_all_sf[k][j] += nhi_vbin
				mf_voronoi_slices_field.phis_all_sf[k][j] += ngal_vbin / mf_voronoi_slices_field.volumes_all_Mpc3[k][j] / mf_voronoi_slices_field.dm



		###############################
		###  repeating for QU galaxies
		###############################

		###  number counting galaxies in massbins
		inds_full = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                   (f.fout.z[f.inds_spatial] >= zlo) &
			                   (f.fout.z[f.inds_spatial] < zhi) &
			                   (f.uvj_class[f.inds_spatial] == 0))
		inds_field = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                    (f.fout.z[f.inds_spatial] >= zlo) &
			                    (f.fout.z[f.inds_spatial] < zhi) &
			                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= 0.3) &
			                    (f.uvj_class[f.inds_spatial] == 0))
		inds_v_lt_0 = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
			                    (f.fout.z[f.inds_spatial] >= zlo) &
			                    (f.fout.z[f.inds_spatial] < zhi) &
			                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= 0) &
			                    (f.uvj_class[f.inds_spatial] == 0))

		digi_mass_full = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_full], mf_voronoi_slices_field.lmassbins)
		ngal_full = pylab.bincount(digi_mass_full, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
		nlo_full, nhi_full = mypy.massfunc.CI(ngal_full)
		mf_voronoi_slices_field.ngals_full_qu[j] += ngal_full
		mf_voronoi_slices_field.elo_ngals_full_qu[j] += nlo_full
		mf_voronoi_slices_field.ehi_ngals_full_qu[j] += nhi_full
		mf_voronoi_slices_field.phis_full_qu[j] += ngal_full / mf_voronoi_slices_field.volumes_full_Mpc3[j] / mf_voronoi_slices_field.dm

		digi_mass_field = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_field], mf_voronoi_slices_field.lmassbins)
		ngal_field = pylab.bincount(digi_mass_field, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
		nlo_field, nhi_field = mypy.massfunc.CI(ngal_field)
		mf_voronoi_slices_field.ngals_field_qu[j] += ngal_field
		mf_voronoi_slices_field.elo_ngals_field_qu[j] += nlo_field
		mf_voronoi_slices_field.ehi_ngals_field_qu[j] += nhi_field
		mf_voronoi_slices_field.phis_field_qu[j] += ngal_field / mf_voronoi_slices_field.volumes_field_Mpc3[j] / mf_voronoi_slices_field.dm

		if len(inds_v_lt_0) > 0:
			digi_mass_v_lt_0 = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_v_lt_0], mf_voronoi_slices_field.lmassbins)
			ngal_v_lt_0 = pylab.bincount(digi_mass_v_lt_0, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
			nlo_v_lt_0, nhi_v_lt_0 = mypy.massfunc.CI(ngal_v_lt_0)
			mf_voronoi_slices_field.ngals_v_lt_0_qu[j] += ngal_v_lt_0
			mf_voronoi_slices_field.elo_ngals_v_lt_0_qu[j] += nlo_v_lt_0
			mf_voronoi_slices_field.ehi_ngals_v_lt_0_qu[j] += nhi_v_lt_0
			mf_voronoi_slices_field.phis_v_lt_0_qu[j] += ngal_v_lt_0 / mf_voronoi_slices_field.volumes_v_lt_0_Mpc3[j] / mf_voronoi_slices_field.dm
		else:
			nlo_v_lt_0, nhi_v_lt_0 = mypy.massfunc.CI(ngal_v_lt_0)
			mf_voronoi_slices_field.elo_ngals_v_lt_0_qu[j] += nlo_v_lt_0
			mf_voronoi_slices_field.ehi_ngals_v_lt_0_qu[j] += nhi_v_lt_0

		###  ... for overdensity bins
		for k in range(len(mf_voronoi_slices_field.overdens_bars)):
			inds_vbin =  pylab.find((f.cat.use[f.inds_spatial] == 1) & 
				                    (f.fout.z[f.inds_spatial] >= zlo) &
				                    (f.fout.z[f.inds_spatial] < zhi) &
				                    (f.uvj_class[f.inds_spatial] == 0) &
				                    (f.voronoi.overdens_matrix[f.inds_spatial, j] >= mf_voronoi_slices_field.overdens_bins[k]) &
				                    (f.voronoi.overdens_matrix[f.inds_spatial, j] <= mf_voronoi_slices_field.overdens_bins[k+1]))

			###  if no galaxies exist in this overdensity bin do nothing
			if len(inds_vbin) == 0:
				pass
			else:
				digi_mass_vbin = pylab.digitize(f.fout.lmass[f.inds_spatial][inds_vbin], mf_voronoi_slices_field.lmassbins)
				ngal_vbin = pylab.bincount(digi_mass_vbin, minlength=len(mf_voronoi_slices_field.lmassbins)+1)[1:-1]
				nlo_vbin, nhi_vbin = mypy.massfunc.CI(ngal_vbin)
				mf_voronoi_slices_field.ngals_all_qu[k][j] += ngal_vbin
				mf_voronoi_slices_field.elo_ngals_all_qu[k][j] += nlo_vbin
				mf_voronoi_slices_field.ehi_ngals_all_qu[k][j] += nhi_vbin
				mf_voronoi_slices_field.phis_all_qu[k][j] += ngal_vbin / mf_voronoi_slices_field.volumes_all_Mpc3[k][j] / mf_voronoi_slices_field.dm


	pickle.dump(mf_voronoi_slices_field, open('../data_after_VoronoiMC_bugfix/MF_%s_voronoi_slices.pickle' % f.name, 'wb'))
	print '\nwrote to: ../data_after_VoronoiMC_bugfix/MF_%s_voronoi_slices.pickle\n' % f.name
















'''

strs = []
strs2 = []
for i in range(len(fields)):
    f = fields[i]
    zlo_phot = f.zclust - 1.5 * f.sigmaz * (1 + f.zclust)
    zhi_phot = f.zclust + 1.5 * f.sigmaz * (1 + f.zclust)
    zlo_spec, zhi_spec = f.zclust, f.zclust
    for ss in f.substructures:
        if ss.zlo_3sig < zlo_spec: zlo_spec = ss.zlo_3sig
        if ss.zhi_3sig > zhi_spec: zhi_spec = ss.zhi_3sig


    members_zspec = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
                               (f.cat.z_spec[f.inds_spatial] > zlo_spec) & 
                               (f.cat.z_spec[f.inds_spatial] < zhi_spec) & 
                               (f.uvj_class[f.inds_spatial] > -1) & 
                               (f.fout.lmass[f.inds_spatial] >= 11.625))

    members_zphot = pylab.find((f.cat.use[f.inds_spatial] == 1) & 
                               (f.cat.z_spec[f.inds_spatial] < 0) & 
                               (f.fout.z[f.inds_spatial] > zlo_phot) & 
                               (f.fout.z[f.inds_spatial] < zhi_phot) & 
                               (f.uvj_class[f.inds_spatial] > -1) & 
                               (f.fout.lmass[f.inds_spatial] >= 11.625))

    str_mugshots = ''
    for i in members_zspec: str_mugshots += ' ./%s/mugshot%05i_%s.pdf' % (f.version, f.cat.id[f.inds_spatial][i], f.version)
    for i in members_zphot: str_mugshots += ' ./%s/mugshot%05i_%s.pdf' % (f.version, f.cat.id[f.inds_spatial][i], f.version)
    strs.append(str_mugshots)

    str_mugshots = ''
    for i in members_zspec: str_mugshots += ' ./%s/mugshot%06i_%s.pdf' % (f.version, f.cat.id[f.inds_spatial][i], f.version)
    for i in members_zphot: str_mugshots += ' ./%s/mugshot%06i_%s.pdf' % (f.version, f.cat.id[f.inds_spatial][i], f.version)
    strs2.append(str_mugshots)


'''















