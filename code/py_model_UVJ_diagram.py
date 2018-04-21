
###  This script is an attempt to test a method of
###  modeling the distribution of galaxies in the
###  UVJ diagram.
###
###  Test 1: Double 2D Gaussian  (FAILED !!!)
###  Test 2: Find trough in the offsets from the diagonal

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
from astropy import wcs, cosmology, constants, units, modeling
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




def uvj_select_region(slope=0.88, intercept=0.59, subplot=False, plot=True, kolor='r', lw=3, ls='-', zorder=99):
    """
    #  DESCRIPTION:
    #    Returns the defined region for selecting
    #    star-forming vs quiescent galaxies from
    #    Whitaker+2011. Will plot it unless told
    #    not to.
    #
    #  INPUT:
    #    z  --->  desired redshift
    #
    #  OUTPUT:
    #    floor  ----->  minimum U-V color to be quiescent
    #    corner1  --->  V-J coordinate at turning point at floor
    #    corner2  --->  U-V coordinate at turning point at wall
    #    wall  ------>  maximum V-J color to be quiescent
    """
    
    floor = 1.3
    wall = 1.6

    c1 = (floor - intercept) / slope
    c2 = slope*wall + intercept

    if plot:
        if subplot:
            ax = subplot.axis()
            p = subplot.plot( (-10,c1), (floor,floor), color=kolor, lw=lw, ls=ls, zorder=zorder )
            p = subplot.plot( (c1,wall), (floor,c2), color=kolor, lw=lw, ls=ls, zorder=zorder )
            p = subplot.plot( (wall,wall), (c2,10), color=kolor, lw=lw, ls=ls, zorder=zorder )
            subplot.axis(ax)
        else:
            ax = pl.axis()
            p = pl.plot( (-10,c1), (floor,floor), color=kolor, lw=lw, ls=ls, zorder=zorder )
            p = pl.plot( (c1,wall), (floor,c2), color=kolor, lw=lw, ls=ls, zorder=zorder )
            p = pl.plot( (wall,wall), (c2,10), color=kolor, lw=lw, ls=ls, zorder=zorder )
            pl.axis(ax)

    return floor, c1, c2, wall




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

data_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'
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
		self.chi2red_thresh = chi2red_thresh

		self.cat = gunzip_read_gzip('%s/%s/%s.cat.gz' % (data_dir, version, version), readcat=1)
		self.zout = gunzip_read_gzip('%s/%s/%s.zout.gz' % (data_dir, version, version), readzout=1)
		self.fout = gunzip_read_gzip('%s/%s/%s.fout.gz' % (data_dir, version, version), readcat=1)
		self.restframe = gunzip_read_gzip('%s/%s/%s.restframe.gz' % (data_dir, version, version), readcat=1)

		self.masslimits = mypy.readcat('../data/completeness/masslimits_%s.dat' % name)

		###  UPDATING USE FLAG WITH REDUCE CHI**2 THRESHOLD
		chi2red = self.zout.chi_p / (self.zout.nfilt - 1.)
		cinds = pylab.find((chi2red > chi2red_thresh) & (self.cat.z_spec < 0))
		self.cat.use[cinds] = 0.


		###  UVJ classification
		self.UV_color = -2.5 * pylab.log10(self.restframe.restflux_U / self.restframe.restflux_V)
		self.VJ_color = -2.5 * pylab.log10(self.restframe.restflux_V / self.restframe.restflux_J)
		self.uvj_class = mypy.uvj_select(self.UV_color, self.VJ_color, self.fout.z)



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
fields.append(field('N200',    'nep200_v0.0.5',       0.691,  0.027, alpha=1./600,  chi2red_thresh=7/3.))
fields.append(field('SC1324',  'sc1324_v0.0.2',       0.755,  0.033, alpha=1./600,  chi2red_thresh=10/3.))
fields.append(field('RCS0224', 'rcs0224-0002_v0.0.2', 0.772,  0.027, alpha=1./500,  chi2red_thresh=10/3.))
fields.append(field('RXJ1716', 'rxj1716+6708_v0.0.8', 0.813,  0.021, alpha=1./500,  chi2red_thresh=8/3.))
fields.append(field('N5281',   'nep5281_v0.0.2',      0.818,  0.029, alpha=1./1000, chi2red_thresh=10/3.))
fields.append(field('SG0023',  'sg0023+0423_v0.1.9',  0.845,  0.025, alpha=1./400,  chi2red_thresh=14/3.))
fields.append(field('SC1604',  'sc1604_v0.0.3',       0.910,  0.029, alpha=1./500,  chi2red_thresh=10/3.))
fields.append(field('SC0910',  'cl0910+5422_v0.0.3',  1.110,  0.035, alpha=1./500,  chi2red_thresh=10/3.))
fields.append(field('SC0849',  'sc0849+4452_v0.0.2',  1.261,  0.029, alpha=1./600,  chi2red_thresh=10/3.))
print ''














outname = '../output/UVJ_offsets.pdf'
pdf = PdfPages(outname)

fig = pylab.figure(figsize=(16., 7))
sp1 = fig.add_subplot(121)
sp2 = fig.add_subplot(122)
fig.subplots_adjust(left=0.08)


for fi in range(len(fields)):

	f = fields[fi]

	sp1.minorticks_on()
	sp2.minorticks_on()
	sp1.axis([-0.5, 2.5, 0., 2.5])

	sp1.set_xlabel('(V - J)$_{\mathrm{rest}}$', size=24)
	sp1.set_ylabel('(U - V)$_{\mathrm{rest}}$', size=24)

	sp2.set_xlabel('offset from diagonal')
	sp2.set_ylabel('Number')

	m_diagonal, b_diagonal = 0.88, 0.59
	m_perp, b_perp = -1. / m_diagonal, 1.6


	###  plotting UVJ points
	zmag_lim = 22.5
	zlo, zhi = 0.6, 1.2
	mypy.uvj_select_region(0.7, sp1, kolor='brown', lw=3)
	inds = pylab.find((f.fout.z > zlo) & 
		              (f.fout.z < zhi) & 
		              (f.zmagnitude < zmag_lim) & 
		              (f.cat.weight_B > 0.05) & 
		              (f.cat.use == 1) &
		              (f.VJ_color > -0.5))

	uv = f.UV_color[inds]
	vj = f.VJ_color[inds]
	xfill = pylab.linspace(-1, 3, 10)
	yfill = m_perp*xfill + b_perp
	sp1.fill_between(xfill, -1, yfill, color='#cccccc')

	asdf0 = sp1.plot(vj, uv, ls='', marker='o', ms=1, mec='r', label='zmag < 22.5\n%.1f < z < %.1f' % (zlo, zhi))
	asdf = sp1.plot([-0.5, 2.5], pylab.array([-0.5, 2.5])*m_diagonal + b_diagonal, color='#666666', lw=3)
	sp1.legend(loc=2, numpoints=4, fontsize=18)



	for p_off in [-0.5, 0.5]:
		y_off = p_off / pylab.sin(pylab.arctan(1. / m_diagonal))
		sp1.plot([-0.5, 2.5], pylab.array([-0.5, 2.5])*m_diagonal + b_diagonal + y_off, color='#666666', ls='--', lw=3)
		sp2.axvline(p_off, color='#666666', ls='--', lw=3)



	###  plotting histogram of offsets from diagonal
	inds = pylab.find((f.fout.z > zlo) & 
		              (f.fout.z < zhi) & 
		              (f.zmagnitude < zmag_lim) & 
		              (f.cat.weight_B > 0.05) & 
		              (f.cat.use == 1) &
		              (f.UV_color > m_perp * f.VJ_color + b_perp))

	uv = f.UV_color[inds]
	vj = f.VJ_color[inds]

	d_uv = (m_diagonal*vj + b_diagonal) - uv
	d_vj = vj - ((uv - b_diagonal) / m_diagonal)
	p = d_vj * d_uv / pylab.sqrt(d_vj**2 + d_uv**2)   # length of perpendicular to diagonal line

	inds_above = pylab.find(uv > m_diagonal*vj+b_diagonal)
	p[inds_above] *= -1

	p_range = (-0.8, 0.8)
	nbins = int(2. * len(inds)**0.5)
	h = sp2.hist(p, histtype='stepfilled', color='#cccccc', lw=1, range=p_range, bins=nbins)


	###  fitting double-Gaussian to offsets hist
	x, y = (h[1][1:] + h[1][:-1]) / 2., h[0]

	double_1D_gaussian = modeling.models.Gaussian1D(y.max(),  0.2, 0.1) + \
	                     modeling.models.Gaussian1D(y.max(), -0.2, 0.1)
	fit_algorithm = modeling.fitting.LevMarLSQFitter()

	gfit = fit_algorithm(double_1D_gaussian, x, y)

	xmod = pylab.linspace(x.min(), x.max(), nbins*10)
	ymod = gfit(xmod)
	sp2.plot(xmod, ymod, color='r')
	sp2.set_xlim(p_range)
	sp2.set_ylim(0, y.max()*1.3)



	###  finding p_offset at the trough of the dounble-gaussian
	x_means = pylab.array([gfit.mean_0.value, gfit.mean_1.value])
	xinds = pylab.find((xmod > x_means.min()) & (xmod < x_means.max()))

	xtrough = xmod[xinds][pylab.find(ymod[xinds] == ymod[xinds].min())[0]]
	sp2.axvline(xtrough, color='r', ls=':', label='offset = %.3f' % xtrough)

	sp2.legend(loc=2, title=f.version)




	pdf.savefig()
	sp1.clear()
	sp2.clear()

pdf.close()
pylab.close()
print '\nwrote to: %s' % outname





