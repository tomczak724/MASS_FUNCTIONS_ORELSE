
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
		self.fout1 = gunzip_read_gzip('%s/%s/%s.fout.gz' % (github_dir, version, version), readcat=1)
		self.fout2 = gunzip_read_gzip('%s/%s/%s_ZFparams.fout.gz' % (github_dir, version, version), readcat=1)


		###  SETTING OBJECTS IDENTIFIED AS SECURE STARS FROM SPECTROSCOPY TO use=0
		self.crossmatch = gunzip_read_gzip('%s/%s/%s.crossmatch.gz' % (github_dir, version, version), readcat=1, dtype=str)
		self.star_inds = pylab.find(self.crossmatch.Q == '-1')
		for i_star in self.star_inds:
			id_phot_arr = self.crossmatch.id_phot[i_star].split(',')
			for id_phot in id_phot_arr:
				if id_phot != '-1':
					self.cat.use[int(id_phot)-1] *= 0


		###  UPDATING USE FLAG WITH REDUCE CHI**2 THRESHOLD
		chi2red = self.zout.chi_p / (self.zout.nfilt - 1.)
		cinds = pylab.find((chi2red > chi2red_thresh) & (self.cat.z_spec < 0))
		self.cat.use[cinds] = 0.


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


dm_fields = []
dm_fields.append(field('N200',    'nep200_v0.0.5',       0.691,  0.027, alpha=1./600, chi2red_thresh=7))
dm_fields.append(field('SC1324',  'sc1324_v0.0.2',       0.755,  0.033, alpha=1./600, chi2red_thresh=10))
dm_fields.append(field('RCS0224', 'rcs0224-0002_v0.0.2', 0.772,  0.027, alpha=1./500, chi2red_thresh=10))
dm_fields.append(field('RXJ1716', 'rxj1716+6708_v0.0.7', 0.813,  0.021, alpha=1./500, chi2red_thresh=8))
dm_fields.append(field('N5281',   'nep5281_v0.0.2',      0.818,  0.029, alpha=1./1000, chi2red_thresh=10))
dm_fields.append(field('SG0023',  'sg0023+0423_v0.1.9',  0.845,  0.025, alpha=1./400, chi2red_thresh=14))
dm_fields.append(field('SC1604',  'sc1604_v0.0.3',       0.910,  0.029, alpha=1./500, chi2red_thresh=10))
dm_fields.append(field('SC0910',  'cl0910+5422_v0.0.3',  1.110,  0.035, alpha=1./500, chi2red_thresh=10))
dm_fields.append(field('SC0849',  'sc0849+4452_v0.0.2',  1.261,  0.029, alpha=1./600, chi2red_thresh=10))
print ''











lmass_b4_2d = []
lmass_afta_2d = []
z_2d = []

lmass_b4_1d = []
lmass_afta_1d = []
z_1d = []

for i in range(len(dm_fields)):

	f = dm_fields[i]

	inds = pylab.find((f.cat.use[f.inds_spatial] == 1) &
		              (f.fout2.lmass[f.inds_spatial] > 8.625) &
		              (f.fout2.z[f.inds_spatial] > 0.45) &
    		          (f.fout2.z[f.inds_spatial] < 1.45))

	lmass_b4_2d.append(f.fout1.lmass[f.inds_spatial][inds])
	lmass_afta_2d.append(f.fout2.lmass[f.inds_spatial][inds])
	z_2d.append(f.fout2.z[f.inds_spatial][inds])

	lmass_b4_1d += f.fout1.lmass[f.inds_spatial][inds].tolist()
	lmass_afta_1d += f.fout2.lmass[f.inds_spatial][inds].tolist()
	z_1d += f.fout2.z[f.inds_spatial][inds].tolist()


lmass_b4_1d = pylab.array(lmass_b4_1d)
lmass_afta_1d = pylab.array(lmass_afta_1d)
z_1d = pylab.array(z_1d)









outname = '../output/deltaMass_ZFparams.pdf'
pdf = PdfPages(outname)


fig = pylab.figure(figsize=(14.3, 6.2))
sp2 = fig.add_subplot(1, 2, 2, aspect='auto')
sp1 = fig.add_subplot(1, 2, 1, aspect='auto')

fig.subplots_adjust(wspace=0, left=0.08)

sp1.set_ylabel('$\Delta$ log(M$_*$)   [ old - new ]', size=22)
sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )', size=22)
sp2.set_xlabel('redshift', size=22)

sp1.grid()
sp2.grid()

sp1.minorticks_on()
sp2.minorticks_on()

ax1, ax2 = [8.6, 11.9, -0.7, 0.7], [0.45, 1.45, -0.7, 0.7]
sp1.axis(ax1)
sp2.axis(ax2)

sp1.axhline(0, color='yellow', lw=4, zorder=2)
sp1.axhline(0, color='k', lw=2, zorder=3)
sp2.axhline(0, color='yellow', lw=4, zorder=2)
sp2.axhline(0, color='k', lw=2, zorder=3)

sp1.axhline(0.25/2, color='yellow', lw=4, ls='-', zorder=2)
sp1.axhline(0.25/2, color='k', lw=2, ls='--', zorder=2)
sp1.axhline(-0.25/2, color='yellow', lw=4, ls='-', zorder=2)
sp1.axhline(-0.25/2, color='k', lw=2, ls='--', zorder=2)

sp2.axhline(0.25/2, color='yellow', lw=4, ls='-', zorder=2)
sp2.axhline(0.25/2, color='k', lw=2, ls='--', zorder=2)
sp2.axhline(-0.25/2, color='yellow', lw=4, ls='-', zorder=2)
sp2.axhline(-0.25/2, color='k', lw=2, ls='--', zorder=2)



#sp1.plot(lmass_afta_1d, lmass_b4_1d - lmass_afta_1d, ls='', marker='o', mec='gray', ms=1, zorder=1)
#sp2.plot(z_1d, lmass_b4_1d - lmass_afta_1d, ls='', marker='o', mec='gray', ms=1, zorder=1)


nbins1, nbins2 = 70, 70
hist2d, xedges, yedges = pylab.histogram2d(lmass_afta_1d, lmass_b4_1d - lmass_afta_1d, bins=(nbins1, nbins2), range=([ax1[0], ax1[1]], [ax1[2], ax1[3]]))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
asdf = sp1.imshow(pylab.log10(hist2d.T+1), extent=extent, interpolation='nearest', cmap=pylab.cm.Greens)
asdf.set_clim(0, pylab.log10(hist2d.max())*1.)
sp1.set_aspect('auto')

nbins1, nbins2 = 70, 70
hist2d, xedges, yedges = pylab.histogram2d(z_1d, lmass_b4_1d - lmass_afta_1d, bins=(nbins1, nbins2), range=([ax2[0], ax2[1]], [ax2[2], ax2[3]]))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
asdf = sp2.imshow(pylab.log10(hist2d.T+1), extent=extent, interpolation='nearest', cmap=pylab.cm.Greens)
asdf.set_clim(0, pylab.log10(hist2d.max())*1.)
sp2.set_aspect('auto')




lmassbins = pylab.arange(8.75, 11.8, 0.25)
lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2
zbins = pylab.arange(0.5, 1.4, 0.1)
zbars = (zbins[1:] + zbins[:-1]) / 2

digi_lmass = pylab.digitize(lmass_afta_1d, lmassbins)
digi_z = pylab.digitize(z_1d, zbins)

dm_medians = []
dm_nmads = []
z_medians = []
z_nmads = []

for dmi in range(1, len(lmassbins)):
	inds = pylab.find(digi_lmass == dmi)
	dm_medians.append(pylab.median(lmass_b4_1d[inds] - lmass_afta_1d[inds]))
	dm_nmads.append(mypy.nmad(lmass_b4_1d[inds] - lmass_afta_1d[inds]))

for zi in range(1, len(zbins)):
	inds = pylab.find(digi_z == zi)
	z_medians.append(pylab.median(lmass_b4_1d[inds] - lmass_afta_1d[inds]))
	z_nmads.append(mypy.nmad(lmass_b4_1d[inds] - lmass_afta_1d[inds]))


sp1.errorbar(lmassbars, dm_medians, yerr=dm_nmads, ls='-', color='r', marker='o', ms=7, mew=2, elinewidth=2, ecolor='r', zorder=5, label='median+nmad')
sp2.errorbar(zbars, z_medians, yerr=z_nmads, ls='-', color='r', marker='o', ms=7, mew=2, elinewidth=2, ecolor='r', zorder=5)

sp1.legend(loc=3, title='All fields', fontsize=16)




pdf.savefig()
sp1.clear()
sp2.clear()








for i in range(len(dm_fields)):

	sp1.set_ylabel('$\Delta$ log(M$_*$)   [ old - new ]', size=22)
	sp1.set_xlabel('log( M$_*$ / M$_{\odot}$ )', size=22)
	sp2.set_xlabel('redshift', size=22)

	sp1.grid()
	sp2.grid()

	sp1.minorticks_on()
	sp2.minorticks_on()

	sp1.axis([8.6, 11.9, -0.7, 0.7])
	sp2.axis([0.45, 1.45, -0.7, 0.7])

	sp1.axhline(0, color='yellow', lw=4, zorder=2)
	sp1.axhline(0, color='k', lw=2, zorder=3)
	sp2.axhline(0, color='yellow', lw=4, zorder=2)
	sp2.axhline(0, color='k', lw=2, zorder=3)

	sp1.axhline(0.25/2, color='yellow', lw=4, ls='-', zorder=2)
	sp1.axhline(0.25/2, color='k', lw=2, ls='--', zorder=2)
	sp1.axhline(-0.25/2, color='yellow', lw=4, ls='-', zorder=2)
	sp1.axhline(-0.25/2, color='k', lw=2, ls='--', zorder=2)

	sp2.axhline(0.25/2, color='yellow', lw=4, ls='-', zorder=2)
	sp2.axhline(0.25/2, color='k', lw=2, ls='--', zorder=2)
	sp2.axhline(-0.25/2, color='yellow', lw=4, ls='-', zorder=2)
	sp2.axhline(-0.25/2, color='k', lw=2, ls='--', zorder=2)



	#sp1.plot(lmass_afta_2d[i], lmass_b4_2d[i] - lmass_afta_2d[i], ls='', marker='o', mec='gray', ms=1, zorder=1)
	#sp2.plot(z_2d[i], lmass_b4_2d[i] - lmass_afta_2d[i], ls='', marker='o', mec='gray', ms=1, zorder=1)

	nbins1, nbins2 = 70, 70
	hist2d, xedges, yedges = pylab.histogram2d(lmass_afta_2d[i], lmass_b4_2d[i] - lmass_afta_2d[i], bins=(nbins1, nbins2), range=([ax1[0], ax1[1]], [ax1[2], ax1[3]]))
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	asdf = sp1.imshow(pylab.log10(hist2d.T+1), extent=extent, interpolation='nearest', cmap=pylab.cm.Greens)
	asdf.set_clim(0, pylab.log10(hist2d.max())*1.)
	sp1.set_aspect('auto')

	nbins1, nbins2 = 70, 70
	hist2d, xedges, yedges = pylab.histogram2d(z_2d[i], lmass_b4_2d[i] - lmass_afta_2d[i], bins=(nbins1, nbins2), range=([ax2[0], ax2[1]], [ax2[2], ax2[3]]))
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	asdf = sp2.imshow(pylab.log10(hist2d.T+1), extent=extent, interpolation='nearest', cmap=pylab.cm.Greens)
	asdf.set_clim(0, pylab.log10(hist2d.max())*1.)
	sp2.set_aspect('auto')



	lmassbins = pylab.arange(8.75, 11.8, 0.25)
	lmassbars = (lmassbins[1:] + lmassbins[:-1]) / 2
	zbins = pylab.arange(0.5, 1.4, 0.1)
	zbars = (zbins[1:] + zbins[:-1]) / 2

	digi_lmass = pylab.digitize(lmass_afta_2d[i], lmassbins)
	digi_z = pylab.digitize(z_2d[i], zbins)

	dm_medians = []
	dm_nmads = []
	z_medians = []
	z_nmads = []

	for dmi in range(1, len(lmassbins)):
		inds = pylab.find(digi_lmass == dmi)
		dm_medians.append(pylab.median(lmass_b4_2d[i][inds] - lmass_afta_2d[i][inds]))
		dm_nmads.append(mypy.nmad(lmass_b4_2d[i][inds] - lmass_afta_2d[i][inds]))

	for zi in range(1, len(zbins)):
		inds = pylab.find(digi_z == zi)
		z_medians.append(pylab.median(lmass_b4_2d[i][inds] - lmass_afta_2d[i][inds]))
		z_nmads.append(mypy.nmad(lmass_b4_2d[i][inds] - lmass_afta_2d[i][inds]))


	sp1.errorbar(lmassbars, dm_medians, yerr=dm_nmads, ls='-', color='r', marker='o', ms=7, mew=2, elinewidth=2, ecolor='r', zorder=5, label='median+nmad')
	sp2.errorbar(zbars, z_medians, yerr=z_nmads, ls='-', color='r', marker='o', ms=7, mew=2, elinewidth=2, ecolor='r', zorder=5)

	sp1.legend(loc=3, title=dm_fields[i].name, fontsize=16)


	pdf.savefig()
	sp1.clear()
	sp2.clear()

pdf.close()
pylab.close()
print '\nwrote to: %s\n' % outname





























###  plotting gridspace b4 and afta

logtau_b4 = pylab.arange(8.5, 10.001, 0.5)
logtau_afta = pylab.arange(7., 11.001, 0.2)

logage_b4 = pylab.arange(8., 10.001, 0.2)
logage_afta = pylab.arange(7.5, 10.1001, 0.1)




pdf = PdfPages('../output/gridspace_FAST.pdf')

fig = pylab.figure(figsize=(13.4, 6.2))

sp1 = fig.add_subplot(121)
sp1.minorticks_on()
sp1.set_xlabel('log( age / yr )')
sp1.set_ylabel(r'log( $\tau$ / yr )')
sp1.axis([7.1, 10.4, 6.5, 11.5])

fig.subplots_adjust(wspace=0, left=0.08)


t = sp1.text(0.05, 0.92, 'Before', transform=sp1.transAxes)

for a in logage_b4:
	sp1.plot([a, a], [logtau_b4[0], logtau_b4[-1]], color='gray')

for t in logtau_b4:
	sp1.plot([logage_b4[0], logage_b4[-1]], [t, t], color='gray')


pdf.savefig()
sp1.clear()








sp1.minorticks_on()
sp1.set_xlabel('log( age / yr )')
sp1.set_ylabel(r'log( $\tau$ / yr )')
sp1.axis([7.1, 10.4, 6.5, 11.5])



t = sp1.text(0.05, 0.92, 'After', transform=sp1.transAxes)

for a in logage_afta:
	sp1.plot([a, a], [logtau_afta[0], logtau_afta[-1]], color='gray')

for t in logtau_afta:
	sp1.plot([logage_afta[0], logage_afta[-1]], [t, t], color='gray')


pdf.savefig()
sp1.clear()
pylab.clf()





sp2 = fig.add_subplot(122)
sp1 = fig.add_subplot(121)
sp1.minorticks_on()
sp2.minorticks_on()
sp1.set_xlabel('log( age / yr )')
sp2.set_xlabel('log( age / yr )')
sp2.set_ylabel(r'log( $\tau$ / yr )')
sp1.set_ylabel(r'log( $\tau$ / yr )')
sp1.axis([7.1, 10.4, 6.5, 11.5])
sp2.axis([7.1, 10.4, 6.5, 11.5])




t = sp1.text(0.05, 0.92, 'Before', transform=sp1.transAxes)
for a in logage_b4:
	sp1.plot([a, a], [logtau_b4[0], logtau_b4[-1]], color='gray')
for t in logtau_b4:
	sp1.plot([logage_b4[0], logage_b4[-1]], [t, t], color='gray')



t = sp2.text(0.05, 0.92, 'After', transform=sp2.transAxes)
for a in logage_afta:
	sp2.plot([a, a], [logtau_afta[0], logtau_afta[-1]], color='gray')
for t in logtau_afta:
	sp2.plot([logage_afta[0], logage_afta[-1]], [t, t], color='gray')




pdf.savefig()
pdf.close()
pylab.close()

print 'wrote to: ../output/gridspace_FAST.pdf'










