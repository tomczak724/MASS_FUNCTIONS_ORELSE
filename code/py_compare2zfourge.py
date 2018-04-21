
import os
import sys
import mypy
import pylab
import subprocess
from scipy import optimize
from mypy import massfunc
from matplotlib.path import Path
from astropy import wcs, cosmology
from threedhst import eazyPy_tomczak
import shapely.geometry as geometry
from shapely.ops import cascaded_union, polygonize
from scipy.spatial import ConvexHull, Delaunay
from matplotlib.backends.backend_pdf import PdfPages

cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

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










data_dir = '/Users/atomczak/GoogleDrive/ORELSE/photometric_catalogs'
data_dir = '/Users/atomczak/GitHub/ORELSE/Catalogs/tomczak_catalogs'

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

class field:
	def __init__(self, name, version, zclust, sigmaz, alpha=50):
		self.name = name          # e.g. "NEP 200"
		self.version = version    # e.g. "nep200_v0.0.4"
		self.zclust = zclust      # cluster redshift
		self.sigmaz = sigmaz      # 1sigma scatter in (zphot-zspec)/(1+zspec)
		self.zspec_lo = 0         # lower redshift bound for specz
		self.zspec_hi = 0         # upper redshift bound for specz
		self.substructures = []

		self.cat = gunzip_read_gzip('%s/%s/%s.cat.gz' % (data_dir, version, version), readcat=1)
		self.zout = gunzip_read_gzip('%s/%s/%s.zout.gz' % (data_dir, version, version), readzout=1)
		self.fout = gunzip_read_gzip('%s/%s/%s.fout.gz' % (data_dir, version, version), readcat=1)


		###  UPDATING USE FLAG WITH REDUCE CHI**2 > 10
		chi2red = self.zout.chi_p / (self.zout.nfilt - 1.)
		cinds = pylab.find(chi2red > 10)
		self.cat.use[cinds] = 0.



		#print '  reading: %s_voronoi.pickle' % name
		#self.voronoi = pickle.load(open('../data/%s_voronoi.pickle' % name, 'rb'))

		xyrd1 = self.cat.x[0], self.cat.y[0], self.cat.ra[0], self.cat.dec[0]
		xyrd2 = self.cat.x[1], self.cat.y[1], self.cat.ra[1], self.cat.dec[1]
		d_arcsec = mypy.radec_sep(xyrd1[2], xyrd1[3], xyrd2[2], xyrd2[3])
		d_pixel = ((xyrd1[0]-xyrd2[0])**2 + (xyrd1[1] - xyrd2[1])**2)**0.5 
		self.px_scale = d_arcsec / d_pixel

		### getting z band magnitude
		try: self.zmagnitude = 25 - 2.5 * pylab.log10(self.cat.fluxauto_z)
		except: pass
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

fields.append(field('N200',    'nep200_v0.0.5',       0.691,  0.027, alpha=1./600))
fields.append(field('SC1324',  'sc1324_v0.0.2',       0.755,  0.033, alpha=1./600))
fields.append(field('RCS0224', 'rcs0224-0002_v0.0.2', 0.772,  0.027, alpha=1./500))
fields.append(field('RXJ1716', 'rxj1716+6708_v0.0.7', 0.813,  0.021, alpha=1./500))
fields.append(field('N5281',   'nep5281_v0.0.2',      0.818,  0.029, alpha=1./1000))
fields.append(field('SG0023',  'sg0023+0423_v0.1.9',  0.845,  0.025, alpha=1./400))
fields.append(field('SC1604',  'sc1604_v0.0.3',       0.910,  0.029, alpha=1./500))
fields.append(field('SC0910',  'cl0910+5422_v0.0.3',  1.110,  0.035, alpha=1./500))
fields.append(field('SC0849',  'sc0849+4452_v0.0.2',  1.261,  0.029, alpha=1./600))



s = ['' for fi in fields]




###  making output directories
try: os.mkdir('../output')
except: pass



dm = 0.25
massbins = pylab.arange(8.-dm/2., 11.5+dm, dm)
massbars = (massbins[1:] + massbins[:-1]) / 2.







###  Field Mass Functions - ZFOURGE 2014
zfourge_dir = '/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables'
zfourge_massfunc_tot = pylab.loadtxt('%s/table1_TOT.dat' % zfourge_dir)
zfourge_massfunc_sf = pylab.loadtxt('%s/table1_SF.dat' % zfourge_dir)
zfourge_massfunc_qui = pylab.loadtxt('%s/table1_QUI.dat' % zfourge_dir)
dm_zfourge = zfourge_massfunc_tot[1][0] - zfourge_massfunc_tot[0][0]



zbins = pylab.array([[0.20, 0.50],
                     [0.50, 0.75],
                     [0.75, 1.00],
                     [1.00, 1.25],
                     [1.25, 1.50],
                     [1.50, 2.00]])
zbars = (zbins[:,1] + zbins[:,0]) / 2.
dz = (zbins[:,1] - zbins[:,0])







pdf = PdfPages('../output/field_mfs/massFunction_comparisons.pdf')

volumes_all = pylab.zeros((len(zbins), len(fields)))
ngal_bins_all = pylab.zeros((len(fields), len(zbins), len(massbars)))





fig = pylab.figure(figsize=(17., 13.5))
sp3 = fig.add_subplot(2, 3, 3)
sp6 = fig.add_subplot(2, 3, 6)
sp2 = fig.add_subplot(2, 3, 2)
sp5 = fig.add_subplot(2, 3, 5)
sp1 = fig.add_subplot(2, 3, 1)
sp4 = fig.add_subplot(2, 3, 4)
sps = [sp1, sp2, sp3, sp4, sp5, sp6]

fig.subplots_adjust(wspace=0, hspace=0, top=0.94, left=0.09, bottom=0.09)






###  plotting Mass Funcitons

for fi in range(len(fields)):

	f = fields[fi]
	title = '%s :   z=%.2f' % (f.version, f.zclust)
	sp2.set_title(title)


	for zi in range(len(zbins)):

		zlo, zhi = zbins[zi]

		volume = massfunc.vcomoving_slice(f.area_arcmin2, zlo, zhi)


		###  integrating EAZY P(z)'s between zphot_lo, zphot_hi
		#integrated_pzs = []
		#zinds = pylab.find((f.zgrid_pz >= zphot_lo) & (f.zgrid_pz <= zphot_hi))
		#for pzi in f.pzs:
			#x = f.zgrid_pz[zinds]
			#y = pzi[zinds]
			#integrated_pzs.append(pylab.trapz(y, x))
		#integrated_pzs = pylab.array(integrated_pzs)






		subinds = pylab.find((f.cat.use[f.inds_spatial] == 1) &
			                 (f.fout.z[f.inds_spatial] > zlo) &
		    	             (f.fout.z[f.inds_spatial] < zhi))

		subinds_massive = pylab.find((f.cat.use[f.inds_spatial] == 1) &
			                         (f.fout.z[f.inds_spatial] > zlo) &
		    	                     (f.fout.z[f.inds_spatial] < zhi) &
		    	                     (f.fout.lmass[f.inds_spatial] > 11))

		if 0.19 < zlo < 0.76:
			for si in subinds_massive:
				s[fi] += ' mugshot_%05i_%s.pdf' % (f.cat.id[f.inds_spatial][si], f.version)







		digi_mass = pylab.digitize(f.fout.lmass[f.inds_spatial][subinds], massbins)

		ngal_bins = pylab.bincount(digi_mass, minlength=len(massbins)+1)[1:-1]
		nlo_poisson, nhi_poisson = [], []
		for n in ngal_bins:
			nhi, nlo = mypy.massfunc.confidence_interval(n)
			nlo_poisson.append(nlo)
			nhi_poisson.append(nhi)
		nlo_poisson, nhi_poisson = pylab.array(nlo_poisson), pylab.array(nhi_poisson)

	
		phi_bins = ngal_bins * 1. / volume / dm

		ephi_lo = phi_bins * nlo_poisson / ngal_bins
		ephi_lo[pylab.isnan(ephi_lo)] = 0
		ephi_hi = phi_bins * nhi_poisson / ngal_bins
		ephi_hi[pylab.isnan(ephi_hi)] = (1. / volume / dm) * nhi_poisson[pylab.isnan(ephi_hi)] / 1.



		###  adding galaxies to total count if the LSS is not in this zbin
		zphot_lo = f.zclust - 1. * f.sigmaz * (1 + f.zclust)
		zphot_hi = f.zclust + 1. * f.sigmaz * (1 + f.zclust)
		if not (zlo < zphot_lo < zhi or zlo < zphot_hi < zhi):
			volumes_all[zi][fi] += volume
			ngal_bins_all[fi][zi] += ngal_bins






		###  plotting it up!!!
		sp = sps[zi]
		sp.grid()
		sp.minorticks_on()
		axis_lims = [7.7, 11.8, -5.4, -0.9]
		sp.axis(axis_lims)
		sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
		sp.set_ylabel('log( $\Phi$ / Mpc$^3$ / dex )', fontsize=22)

		if zlo < f.zclust < zhi:
			sp.fill_between(axis_lims[0:2], axis_lims[2], axis_lims[3], color='#ffe6e6', zorder=0)


		sp.text(0.5, 0.9, '%.2f < z < %.2f' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')


		###  plotting zfourge first
		i_mf = 3*zi + 1
		inds = pylab.find(zfourge_massfunc_tot[:, i_mf] > -98)
		sp.fill_between(zfourge_massfunc_tot[inds, 0],
			            zfourge_massfunc_tot[inds, i_mf] + zfourge_massfunc_tot[inds, i_mf+1], 
			            zfourge_massfunc_tot[inds, i_mf] - zfourge_massfunc_tot[inds, i_mf+2], 
			            color='#bfbfbf', zorder=1)
		sp.axvline(0, color='#bfbfbf', lw=10, label='ZFOURGE 2014')
	

		###  plotting orelse next
		y = pylab.log10(phi_bins)
		ylo = pylab.log10(phi_bins) - pylab.log10(phi_bins - ephi_lo)
		yhi = pylab.log10(phi_bins + ephi_hi) - pylab.log10(phi_bins)
		sp.errorbar(massbars, y, xerr=dm/2., yerr=[ylo, yhi], ls='', marker='o', ms=7, mew=1.5, \
			        elinewidth=1.5, mfc='#cc00ff', ecolor='k', label='ORELSE today', zorder=2)


		###  fitting double Schechter
		#yerr = (ylo + yhi) / 2.
		#fit, cov = optimize.curve_fit(massfunc.dschechter_mf, massbars, 10**y, p0=[10.78, -0.98, -2.54, -1.9, -4.29], sigma=yerr)

		#xmodel = pylab.linspace(9, 11.7, 1000)
		#ymodel = pylab.log10(massfunc.dschechter_mf(xmodel, *fit))
		#sp.plot(xmodel, ymodel, color='r', label='2x Schechter fit', zorder=1)


		if zi == 3:
			sp.legend(loc=3, fontsize=18, title=title)


	pylab.savefig('../output/field_mfs/massFunction_comparisons_%s.png' % f.name)
	pdf.savefig()
	for sp in sps:
		sp.clear()



###  plotting all fields averaged together
dschechter_fits = []
dschechter_errs = []
for zi in range(len(zbins)):

	zlo, zhi = zbins[zi]
	title = 'All fields minus LSS'



	ngal = ngal_bins_all.sum(axis=0)[zi]
	vol = volumes_all[zi].sum()
	nlo_poisson, nhi_poisson = [], []
	for n in ngal:
		nhi, nlo = mypy.massfunc.confidence_interval(n)
		nlo_poisson.append(nlo)
		nhi_poisson.append(nhi)
	nlo_poisson, nhi_poisson = pylab.array(nlo_poisson), pylab.array(nhi_poisson)

	phi = ngal * 1. / vol / dm

	ephi_lo = phi * nlo_poisson / ngal
	ephi_lo[pylab.isnan(ephi_lo)] = 0
	ephi_hi = phi * nhi_poisson / ngal
	ephi_hi[pylab.isnan(ephi_hi)] = (1. / vol / dm) * nhi_poisson[pylab.isnan(ephi_hi)] / 1.





	sp = sps[zi]
	sp.grid()
	sp.minorticks_on()
	axis_lims = [7.7, 11.8, -5.4, -0.9]
	sp.axis(axis_lims)
	sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
	sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
	sp.set_ylabel('log( $\Phi$ / Mpc$^3$ / dex )', fontsize=22)


	sp.text(0.5, 0.9, '%.2f < z < %.2f' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')


	###  plotting zfourge first
	i_mf = 3*zi + 1
	inds = pylab.find(zfourge_massfunc_tot[:, i_mf] > -98)
	sp.fill_between(zfourge_massfunc_tot[inds, 0],
		            zfourge_massfunc_tot[inds, i_mf] + zfourge_massfunc_tot[inds, i_mf+1], 
		            zfourge_massfunc_tot[inds, i_mf] - zfourge_massfunc_tot[inds, i_mf+2], 
		            color='#bfbfbf', zorder=1)
	sp.axvline(0, color='#bfbfbf', lw=10, label='ZFOURGE 2014')


	###  plotting orelse next
	y = pylab.log10(phi)
	ylo = pylab.log10(phi) - pylab.log10(phi - ephi_lo)
	yhi = pylab.log10(phi + ephi_hi) - pylab.log10(phi)
	sp.errorbar(massbars, y, xerr=dm/2., yerr=[ylo, yhi], ls='', marker='o', ms=7, mew=1.5, \
		        elinewidth=1.5, mfc='k', ecolor='k', label='ORELSE today', zorder=2)


	###  fitting double Schechter
	mlim = 9.25
	inds_fit = pylab.find(massbars >= mlim)
	yerr = (ylo + yhi) / 2.
	fit, cov = optimize.curve_fit(massfunc.dschechter_mf, massbars[inds_fit], 10**y[inds_fit], p0=[10.78, -0.98, -2.54, -1.9, -4.29], sigma=(ephi_hi[inds_fit]+ephi_lo[inds_fit]))
	dschechter_fits.append(fit)
	dschechter_errs.append(cov.diagonal()**0.5)

	xmodel = pylab.linspace(mlim, 11.7, 1000)
	ymodel = pylab.log10(massfunc.dschechter_mf(xmodel, *fit))
	sp.plot(xmodel, ymodel, color='r', label='2x Schechter fit', zorder=1)


	if zi == 3:
		sp.legend(loc=3, fontsize=18, title=title)


pylab.savefig('../output/field_mfs/massFunction_comparisons_combine.png')
pdf.savefig()


pdf.close()
pylab.close()

























###  plotting schechter params vs. redshift

dschechter_fits = pylab.array(dschechter_fits)
dschechter_errs = pylab.array(dschechter_errs)

zfourge_params = mypy.readcat('/Users/atomczak/PROJECTS/MASS-FUNCTIONS_v2/tables/table2_TOT.dat')
inds = pylab.arange(len(zbins))


fig = pylab.figure(figsize=(23, 8.5))
sp1 = fig.add_subplot(131)
sp2 = fig.add_subplot(132)
sp3 = fig.add_subplot(133)

sp1.minorticks_on()
sp2.minorticks_on()
sp3.minorticks_on()

sp1.grid()
sp2.grid()
sp3.grid()

fig.subplots_adjust(top=0.94, left=0.08)



###  Mstar
sp1.set_xlabel('redshfit')
sp1.set_ylabel('log( M$^*$ / M$_{\odot}$ )')
sp1.axis([-0.2, 2.3, 10.35, 11.4])

y = dschechter_fits[:,0]
yerr = dschechter_errs[:,0]
sp1.errorbar(zbars + 0.015, y, xerr=dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
	         ecolor='r', elinewidth=2, label='ORELSE today')


sp1.errorbar(zbars - 0.015, zfourge_params.Mstar[inds], xerr=dz/2., yerr=zfourge_params.e_Mstar[inds], \
	         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')






###  alpha2
sp2.set_xlabel('redshfit')
sp2.set_ylabel('alpha 2')
sp2.axis([-0.2, 2.3, -2.2, -1])

y = dschechter_fits[:,3]
yerr = dschechter_errs[:,3]
sp2.errorbar(zbars + 0.015, y, xerr=dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
	         ecolor='r', elinewidth=2, label='ORELSE today')


sp2.errorbar(zbars - 0.015, zfourge_params.alpha2[inds], xerr=dz/2., yerr=zfourge_params.e_alpha2[inds], \
	         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')







###  log(phistar1 + phistar2)
sp3.set_xlabel('redshfit')
sp3.set_ylabel('log( $\Phi_1^*$ + $\Phi_2^*$ )')
sp3.axis([-0.2, 2.3, -3.4, -2.1])

y = pylab.log10(10**dschechter_fits[:,2] + 10**dschechter_fits[:,4])

dphi1 = (10**dschechter_errs[:,2] - 1) * 10**dschechter_fits[:,2]
dphi2 = (10**dschechter_errs[:,4] - 1) * 10**dschechter_fits[:,4]
yerr = pylab.sqrt((dphi1**2 + dphi2**2) / pylab.log(10)**2 / (10**dschechter_fits[:,2] + 10**dschechter_fits[:,4])**2)

sp3.errorbar(zbars + 0.015, y, xerr=dz/2., yerr=yerr, ls='', marker='o', ms=9, mfc='r', \
	         ecolor='r', elinewidth=2, label='ORELSE today')


y = pylab.log10(10**zfourge_params.phistar1 + 10**zfourge_params.phistar2)

dphi1 = (10**zfourge_params.e_phistar1 - 1) * 10**zfourge_params.phistar1
dphi2 = (10**zfourge_params.e_phistar2 - 1) * 10**zfourge_params.phistar2
yerr = pylab.sqrt((dphi1**2 + dphi2**2) / pylab.log(10)**2 / (10**zfourge_params.phistar1 + 10**zfourge_params.phistar2)**2)

sp3.errorbar(zbars - 0.015, y[inds], xerr=dz/2., yerr=yerr[inds], \
	         ls='', marker='s', ms=9, mfc='k', ecolor='k', elinewidth=2, label='ZFOURGE 2014')

sp3.legend(loc=3, fontsize=18)


fig.savefig('../output/field_mfs/schechter_params.png')
pylab.close()






















###  plotting MF Residuals wrt ZFOURGE
pdf = PdfPages('../output/field_mfs/massFunction_residuals.pdf')

fig = pylab.figure(figsize=(17., 13.5))
sp3 = fig.add_subplot(2, 3, 3)
sp6 = fig.add_subplot(2, 3, 6)
sp2 = fig.add_subplot(2, 3, 2)
sp5 = fig.add_subplot(2, 3, 5)
sp1 = fig.add_subplot(2, 3, 1)
sp4 = fig.add_subplot(2, 3, 4)
sps = [sp1, sp2, sp3, sp4, sp5, sp6]

fig.subplots_adjust(wspace=0, hspace=0, top=0.94, left=0.09, bottom=0.09)






###  plotting Mass Funcitons

for fi in range(len(fields)):

	f = fields[fi]
	title = '%s :   z=%.2f' % (f.version, f.zclust)
	sp2.set_title(title)


	for zi in range(len(zbins)):

		zlo, zhi = zbins[zi]

		volume = mypy.massfunc.vcomoving_slice(f.area_arcmin2, zlo, zhi)


		###  integrating EAZY P(z)'s between zphot_lo, zphot_hi
		#integrated_pzs = []
		#zinds = pylab.find((f.zgrid_pz >= zphot_lo) & (f.zgrid_pz <= zphot_hi))
		#for pzi in f.pzs:
			#x = f.zgrid_pz[zinds]
			#y = pzi[zinds]
			#integrated_pzs.append(pylab.trapz(y, x))
		#integrated_pzs = pylab.array(integrated_pzs)






		subinds = pylab.find((f.cat.use[f.inds_spatial] == 1) &
			                 (f.fout.z[f.inds_spatial] > zlo) &
		    	             (f.fout.z[f.inds_spatial] < zhi))







		digi_mass = pylab.digitize(f.fout.lmass[f.inds_spatial][subinds], massbins)

		ngal_bins = pylab.bincount(digi_mass, minlength=len(massbins)+1)[1:-1]
		nlo_poisson, nhi_poisson = [], []
		for n in ngal_bins:
			nhi, nlo = mypy.massfunc.confidence_interval(n)
			nlo_poisson.append(nlo)
			nhi_poisson.append(nhi)
		nlo_poisson, nhi_poisson = pylab.array(nlo_poisson), pylab.array(nhi_poisson)

	
		phi_bins = ngal_bins * 1. / volume / dm

		ephi_lo = phi_bins * nlo_poisson / ngal_bins
		ephi_lo[pylab.isnan(ephi_lo)] = 0
		ephi_hi = phi_bins * nhi_poisson / ngal_bins
		ephi_hi[pylab.isnan(ephi_hi)] = (1. / volume / dm) * nhi_poisson[pylab.isnan(ephi_hi)] / 1.








		###  plotting it up!!!
		sp = sps[zi]
		sp.grid()
		sp.minorticks_on()
		axis_lims = [7.7, 11.8, -0.7, 1.4]
		sp.axis(axis_lims)
		sp.axhline(0, color='k', lw=1)
		sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
		sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
		sp.set_ylabel('log( $\Phi_{orelse}$ / $\Phi_{zfourge}$ )', fontsize=22)

		if zlo < f.zclust < zhi:
			sp.fill_between(axis_lims[0:2], axis_lims[2], axis_lims[3], color='#ffe6e6', zorder=0)


		sp.text(0.5, 0.9, '%.2f < z < %.2f' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')



		###  plotting zfourge first
		i_mf = 3*zi + 1
		inds = pylab.find(zfourge_massfunc_tot[:, i_mf] > -98)
		sp.fill_between(zfourge_massfunc_tot[inds, 0],
			            zfourge_massfunc_tot[inds, i_mf+1], 
			            -zfourge_massfunc_tot[inds, i_mf+2], 
			            color='#bfbfbf', zorder=1)
		sp.axvline(0, color='#bfbfbf', lw=10, label='ZFOURGE 2014')
	

		###  plotting orelse next
		y = pylab.log10(phi_bins[inds]) - zfourge_massfunc_tot[inds, i_mf]
		ylo = pylab.log10(phi_bins[inds]) - pylab.log10(phi_bins[inds] - ephi_lo[inds])
		yhi = pylab.log10(phi_bins[inds] + ephi_hi[inds]) - pylab.log10(phi_bins[inds])
		sp.errorbar(massbars[inds], y, xerr=dm/2., yerr=[ylo, yhi], ls='', marker='o', ms=7, mew=1.5, \
			        elinewidth=1.5, mfc='#cc00ff', ecolor='k', label='ORELSE today', zorder=2)

		if zi == 3:
			leg = sp.legend(loc=3, fontsize=15, title=title)
			leg.get_title().set_fontsize(15)

	pylab.savefig('../output/field_mfs/massFunction_residuals_%s.png' % f.name)
	pdf.savefig()
	for sp in sps:
		sp.clear()



###  plotting all fields averaged together
for zi in range(len(zbins)):

	zlo, zhi = zbins[zi]
	title = 'All fields minus LSS'



	ngal = ngal_bins_all.sum(axis=0)[zi]
	vol = volumes_all[zi].sum()
	nlo_poisson, nhi_poisson = [], []
	for n in ngal:
		nhi, nlo = mypy.massfunc.confidence_interval(n)
		nlo_poisson.append(nlo)
		nhi_poisson.append(nhi)
	nlo_poisson, nhi_poisson = pylab.array(nlo_poisson), pylab.array(nhi_poisson)

	phi = ngal * 1. / vol / dm

	ephi_lo = phi * nlo_poisson / ngal
	ephi_lo[pylab.isnan(ephi_lo)] = 0
	ephi_hi = phi * nhi_poisson / ngal
	ephi_hi[pylab.isnan(ephi_hi)] = (1. / vol / dm) * nhi_poisson[pylab.isnan(ephi_hi)] / 1.





	sp = sps[zi]
	sp.grid()
	sp.minorticks_on()
	axis_lims = [7.7, 11.8, -0.7, 1.4]
	sp.axis(axis_lims)
	sp.axhline(0, color='k', lw=1)
	sp.set_xticklabels(['', '8', '', '9', '', '10', '', '11'])
	sp.set_xlabel('log( M$_*$ / M$_{\odot}$ )', fontsize=22)
	sp.set_ylabel('log( $\Phi_{orelse}$ / $\Phi_{zfourge}$ )', fontsize=22)


	sp.text(0.5, 0.9, '%.2f < z < %.2f' % (zlo, zhi), transform=sp.transAxes, fontweight='bold')


	###  plotting zfourge first
	i_mf = 3*zi + 1
	inds = pylab.find(zfourge_massfunc_tot[:, i_mf] > -98)
	sp.fill_between(zfourge_massfunc_tot[inds, 0],
		            zfourge_massfunc_tot[inds, i_mf+1], 
		            -zfourge_massfunc_tot[inds, i_mf+2], 
		            color='#bfbfbf', zorder=1)
	sp.axvline(0, color='#bfbfbf', lw=10, label='ZFOURGE 2014')


	###  plotting orelse next
	y = pylab.log10(phi[inds]) - zfourge_massfunc_tot[inds, i_mf]
	ylo = pylab.log10(phi[inds]) - pylab.log10(phi[inds] - ephi_lo[inds])
	yhi = pylab.log10(phi[inds] + ephi_hi[inds]) - pylab.log10(phi[inds])
	sp.errorbar(massbars[inds], y, xerr=dm/2., yerr=[ylo, yhi], ls='', marker='o', ms=7, mew=1.5, \
		        elinewidth=1.5, mfc='k', ecolor='k', label='ORELSE today', zorder=2)

	if zi == 3:
		leg = sp.legend(loc=3, fontsize=15, title=title)
		leg.get_title().set_fontsize(15)

pylab.savefig('../output/field_mfs/massFunction_residuals_combine.png')
pdf.savefig()


pdf.close()
pylab.close()






















